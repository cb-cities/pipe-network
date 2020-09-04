#include "matrix_assembler.h"

pipenetwork::linear_system::MatrixAssembler::MatrixAssembler(
    const std::shared_ptr<pipenetwork::Mesh>& mesh,
    std::shared_ptr<pipenetwork::Curves>& curves_info, bool pdd_mode)
    : mesh_{mesh}, curves_info_{curves_info}, pdd_{pdd_mode} {
  vars_ = std::make_shared<Variables>(mesh_);
  update_curves_info_();
  res_ = std::make_shared<Residuals>(mesh_, vars_, curves_info_);
  jac_ = std::make_shared<Jacobian>(mesh_, vars_, res_, curves_info_);
}

void pipenetwork::linear_system::MatrixAssembler::assemble_residual() {
  if (pdd_) {
    res_->assemble_residual_pdd();
  } else {
    res_->assemble_residual();
  }
}

void pipenetwork::linear_system::MatrixAssembler::update_jacobian() {
  if (pdd_) {
    jac_->update_jacobian_pdd();
  } else {
    jac_->update_jacobian();
  }
}

void pipenetwork::linear_system::MatrixAssembler::update_curves_info_() {
  auto leak_ids = mesh_->leak_nids();
  auto leak_areas = vars_->leak_areas();
  for (Index i = 0; i < leak_ids.size(); i++) {
    std::string leak_id = std::to_string(leak_ids[i]);
    auto leak_area = leak_areas[i];
    curves_info_->add_leak_poly_vec(leak_id, leak_area);
  }
}
