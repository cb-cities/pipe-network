#include "cuda_unsym.h"
#include "factory.h"

using namespace std;
static Register<pipenetwork::Solver, pipenetwork::Cuda_unsym> registry("cuda");

Eigen::VectorXd pipenetwork::Cuda_unsym::solve() {
  double start, stop;
  double time_solve_gpu;
  assign_mtx_ptrs();
  allocate_host_memories();
  allocate_device_memories();
  init_cuda_handlers();

  //    reorder_mtx_A();
  //    create_mtx_B();
  //    create_mtx_QB();
  cp_hdata_ddata();
  //    solve_z();
  solve_x();
  x_to_host();
  create_x_diff();
  release_resources();

  return x_diff_;
}

void pipenetwork::Cuda_unsym::assign_mtx_ptrs() {
  h_csrValA_ = a_;
  h_csrRowPtrA_ = ia_;
  h_csrColIndA_ = ja_;
  baseA_ = h_csrRowPtrA_[0];
  h_b_ = vec_b_->data();
}

void pipenetwork::Cuda_unsym::allocate_host_memories() {
  h_z_ = (double*)malloc(sizeof(double) * colsA_);
  h_x_ = (double*)malloc(sizeof(double) * colsA_);
  h_Qb_ = (double*)malloc(sizeof(double) * rowsA_);
  h_r_ = (double*)malloc(sizeof(double) * rowsA_);

  h_Q_ = (int*)malloc(sizeof(int) * colsA_);
  h_csrRowPtrB_ = (int*)malloc(sizeof(int) * (rowsA_ + 1));
  h_csrColIndB_ = (int*)malloc(sizeof(int) * nnzA_);
  h_csrValB_ = (double*)malloc(sizeof(double) * nnzA_);
  h_mapBfromA_ = (int*)malloc(sizeof(int) * nnzA_);

  assert(NULL != h_z_);
  assert(NULL != h_x_);
  assert(NULL != h_b_);
  assert(NULL != h_Qb_);
  assert(NULL != h_r_);
  assert(NULL != h_Q_);
  assert(NULL != h_csrRowPtrB_);
  assert(NULL != h_csrColIndB_);
  assert(NULL != h_csrValB_);
  assert(NULL != h_mapBfromA_);
}

void pipenetwork::Cuda_unsym::allocate_device_memories() {
  checkCudaErrors(
      cudaMalloc((void**)&d_csrRowPtrA_, sizeof(int) * (rowsA_ + 1)));
  checkCudaErrors(cudaMalloc((void**)&d_csrColIndA_, sizeof(int) * nnzA_));
  checkCudaErrors(cudaMalloc((void**)&d_csrValA_, sizeof(double) * nnzA_));
  checkCudaErrors(
      cudaMalloc((void**)&d_csrRowPtrB_, sizeof(int) * (rowsA_ + 1)));
  checkCudaErrors(cudaMalloc((void**)&d_csrColIndB_, sizeof(int) * nnzA_));
  checkCudaErrors(cudaMalloc((void**)&d_csrValB_, sizeof(double) * nnzA_));
  checkCudaErrors(cudaMalloc((void**)&d_Q_, sizeof(int) * colsA_));
  checkCudaErrors(cudaMalloc((void**)&d_z_, sizeof(double) * colsA_));
  checkCudaErrors(cudaMalloc((void**)&d_x_, sizeof(double) * colsA_));
  checkCudaErrors(cudaMalloc((void**)&d_b_, sizeof(double) * rowsA_));
  checkCudaErrors(cudaMalloc((void**)&d_Qb_, sizeof(double) * rowsA_));
  checkCudaErrors(cudaMalloc((void**)&d_r_, sizeof(double) * rowsA_));
}

void pipenetwork::Cuda_unsym::init_cuda_handlers() {
  checkCudaErrors(cusolverSpCreate(&handle_));
  checkCudaErrors(cusparseCreate(&cusparseHandle_));

  checkCudaErrors(cudaStreamCreate(&stream_));
  /* bind stream to cusparse and cusolver*/
  checkCudaErrors(cusolverSpSetStream(handle_, stream_));
  checkCudaErrors(cusparseSetStream(cusparseHandle_, stream_));

  /* configure matrix descriptor*/
  checkCudaErrors(cusparseCreateMatDescr(&descrA_));
  checkCudaErrors(cusparseSetMatType(descrA_, CUSPARSE_MATRIX_TYPE_GENERAL));
  if (baseA_) {
    checkCudaErrors(cusparseSetMatIndexBase(descrA_, CUSPARSE_INDEX_BASE_ONE));
  } else {
    checkCudaErrors(cusparseSetMatIndexBase(descrA_, CUSPARSE_INDEX_BASE_ZERO));
  }
}

void pipenetwork::Cuda_unsym::reorder_mtx_A() {

  checkCudaErrors(cusolverSpXcsrmetisndHost(handle_, rowsA_, nnzA_, descrA_,
                                            h_csrRowPtrA_, h_csrColIndA_,
                                            NULL, /* default setting. */
                                            h_Q_));
}

void pipenetwork::Cuda_unsym::create_mtx_B() {
  // step 2.2: B = A(Q,Q)
  memcpy(h_csrRowPtrB_, h_csrRowPtrA_, sizeof(int) * (rowsA_ + 1));
  memcpy(h_csrColIndB_, h_csrColIndA_, sizeof(int) * nnzA_);

  checkCudaErrors(cusolverSpXcsrperm_bufferSizeHost(
      handle_, rowsA_, colsA_, nnzA_, descrA_, h_csrRowPtrB_, h_csrColIndB_,
      h_Q_, h_Q_, &size_perm_));
  if (buffer_cpu_) {
    free(buffer_cpu_);
  }
  buffer_cpu_ = (void*)malloc(sizeof(char) * size_perm_);
  assert(NULL != buffer_cpu);

  /* h_mapBfromA = Identity */
  for (int j = 0; j < nnzA_; j++) {
    h_mapBfromA_[j] = j;
  }
  checkCudaErrors(cusolverSpXcsrpermHost(
      handle_, rowsA_, colsA_, nnzA_, descrA_, h_csrRowPtrB_, h_csrColIndB_,
      h_Q_, h_Q_, h_mapBfromA_, buffer_cpu_));

  /* B = A( mapBfromA ) */
  for (int j = 0; j < nnzA_; j++) {
    h_csrValB_[j] = h_csrValA_[h_mapBfromA_[j]];
  }
}

void pipenetwork::Cuda_unsym::create_mtx_QB() {
  /* h_Qb = b(Q) */
  for (int row = 0; row < rowsA_; row++) {
    h_Qb_[row] = h_b_[h_Q_[row]];
  }
}

void pipenetwork::Cuda_unsym::cp_hdata_ddata() {
  checkCudaErrors(cudaMemcpyAsync(d_csrRowPtrA_, h_csrRowPtrA_,
                                  sizeof(int) * (rowsA_ + 1),
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_csrColIndA_, h_csrColIndA_,
                                  sizeof(int) * nnzA_, cudaMemcpyHostToDevice,
                                  stream_));
  checkCudaErrors(cudaMemcpyAsync(d_csrValA_, h_csrValA_,
                                  sizeof(double) * nnzA_,
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_csrRowPtrB_, h_csrRowPtrB_,
                                  sizeof(int) * (rowsA_ + 1),
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_csrColIndB_, h_csrColIndB_,
                                  sizeof(int) * nnzA_, cudaMemcpyHostToDevice,
                                  stream_));
  checkCudaErrors(cudaMemcpyAsync(d_csrValB_, h_csrValB_,
                                  sizeof(double) * nnzA_,
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_b_, h_b_, sizeof(double) * rowsA_,
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_Qb_, h_Qb_, sizeof(double) * rowsA_,
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaMemcpyAsync(d_Q_, h_Q_, sizeof(int) * rowsA_,
                                  cudaMemcpyHostToDevice, stream_));
  checkCudaErrors(cudaDeviceSynchronize());
}

void pipenetwork::Cuda_unsym::solve_z() {
  checkCudaErrors(cusolverSpDcsrlsvqr(
      handle_, rowsA_, nnzA_, descrA_, d_csrValB_, d_csrRowPtrB_, d_csrColIndB_,
      d_Qb_, tol_, reorder_, d_z_, &singularity_));

  checkCudaErrors(cudaDeviceSynchronize());
  if (0 <= singularity_) {
    printf("WARNING: the matrix is singular at row %d under tol (%E)\n",
           singularity_, tol_);
  }
};

void pipenetwork::Cuda_unsym::solve_x() {

  checkCudaErrors(cusolverSpDcsrlsvqr(
      handle_, rowsA_, nnzA_, descrA_, d_csrValA_, d_csrRowPtrA_, d_csrColIndA_,
      d_b_, tol_, reorder_, d_x_, &singularity_));

  checkCudaErrors(cudaDeviceSynchronize());
  if (0 <= singularity_) {
    printf("WARNING: the matrix is singular at row %d under tol (%E)\n",
           singularity_, tol_);
  }

  //    checkCudaErrors(cusparseDsctr(cusparseHandle_,
  //                                  rowsA_,
  //                                  d_z_,
  //                                  d_Q_,
  //                                  d_x_,
  //                                  CUSPARSE_INDEX_BASE_ZERO));
  //    checkCudaErrors(cudaDeviceSynchronize());
}
void pipenetwork::Cuda_unsym::x_to_host() {
  checkCudaErrors(cudaMemcpyAsync(h_x_, d_x_, sizeof(double) * colsA_,
                                  cudaMemcpyDeviceToHost, stream_));
  checkCudaErrors(cudaDeviceSynchronize());
}
void pipenetwork::Cuda_unsym::create_x_diff() {
  x_diff_.resize(colsA_);
  for (int i = 0; i < colsA_; i++) {
    x_diff_[i] = h_x_[i];
  }
}
void pipenetwork::Cuda_unsym::release_resources() {
  release_host_resources();
  release_device_resources();
}
void pipenetwork::Cuda_unsym::release_host_resources() {
  if (handle_) {
    checkCudaErrors(cusolverSpDestroy(handle_));
  }
  if (cusparseHandle_) {
    checkCudaErrors(cusparseDestroy(cusparseHandle_));
  }
  if (stream_) {
    checkCudaErrors(cudaStreamDestroy(stream_));
  }
  if (descrA_) {
    checkCudaErrors(cusparseDestroyMatDescr(descrA_));
  }

  if (h_z_) {
    free(h_z_);
  }
  if (h_x_) {
    free(h_x_);
  }
  if (h_Qb_) {
    free(h_Qb_);
  }
  if (h_r_) {
    free(h_r_);
  }

  if (h_Q_) {
    free(h_Q_);
  }

  //    if (buffer_cpu_) { free(buffer_cpu_); }
}
void pipenetwork::Cuda_unsym::release_device_resources() {
  if (d_csrValA_) {
    checkCudaErrors(cudaFree(d_csrValA_));
  }
  if (d_csrRowPtrA_) {
    checkCudaErrors(cudaFree(d_csrRowPtrA_));
  }
  if (d_csrColIndA_) {
    checkCudaErrors(cudaFree(d_csrColIndA_));
  }
  if (d_csrValB_) {
    checkCudaErrors(cudaFree(d_csrValB_));
  }
  if (d_csrRowPtrB_) {
    checkCudaErrors(cudaFree(d_csrRowPtrB_));
  }
  if (d_csrColIndB_) {
    checkCudaErrors(cudaFree(d_csrColIndB_));
  }
  if (d_Q_) {
    checkCudaErrors(cudaFree(d_Q_));
  }
  if (d_z_) {
    checkCudaErrors(cudaFree(d_z_));
  }
  if (d_x_) {
    checkCudaErrors(cudaFree(d_x_));
  }
  if (d_b_) {
    checkCudaErrors(cudaFree(d_b_));
  }
  if (d_Qb_) {
    checkCudaErrors(cudaFree(d_Qb_));
  }
  if (d_r_) {
    checkCudaErrors(cudaFree(d_r_));
  }
}
