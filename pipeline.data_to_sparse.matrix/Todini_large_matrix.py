#!Python3
### Based on the Global Algorithm Representation in Todini and Rossman (2013) (Eq 16 in the paper).
### According to the paper, for a water pipeline network with N nodes and p edges, no flow is observed. All pipe flows are unknown (the unknown flow's pipe index number stored in list varaible uQ). At least one nodal head is known. The index of the unknown heads are stored in the list variable uH. The number of unknown heads (nuH) satisfies 0 < nuH < N.
### The number of equations (mass conservation and energy conservation) equal to the number of unknowns (nuQ + nuH). Thus the system of equations can be solved using scipy.optimize.fsolve.

import numpy as np
from scipy.optimize import fsolve
import scipy.sparse
import sys

### Test with small network
# r = [0.000298867,0.007789731,0.014060475,0.002157818, 0.011684596,0.02812095,0.014060475]
# Ab12 = np.array([                  #fixed topology, p x N matrix
#         [-1, 1, 0, 0, 0],
#         [-1, 0, 1, 0, 0],
#         [ 0,-1, 1, 0, 0],
#         [ 0,-1, 0, 1, 0],
#         [ 0,-1, 0, 0, 1],
#         [ 0, 0,-1, 0, 1],
#         [ 0, 0, 0,-1, 1]])
# Ab12 = scipy.sparse.csr_matrix(Ab12)
# q_of_uH = [-10, -20, -30, -40]
### End

a = 1.852       # the alpha coefficient Hazen-Williams
r = 0.0003      # the r_k coefficient that depends on the dimensions used as well as on the pipe diameter, length, and possibly the flow rate

Ab12 = scipy.sparse.load_npz('sparse_matrix_bz247.npz')   # A_bar_12, the topology sparse matrix, p*N
(p, N) = Ab12.shape
print('# edges:', p, '# nodes', N)

### Now, suppose only one hydraulic head, the one of node 0, is known (H0=100)
### The demand from node 1 is -100. Demands for node 2 to N are 0. (the q_of_uH list)
### All flows are unknown. The demand of node 0 as well as the hydraulic head for node 1 to N are unknown.
                                     
A10 = Ab12[:,0]                      # creation of matrix A10, the first column of Ab12, corresponding the topology relationship for the Node 0.
H0 = [100]                 # list of known head, here the head of Node 0 is 100.

nuH = N-1                            # number of unknown heads. In this case the head from node 1 to N are all unknown.
A12 = Ab12[:,1:N]                    # creation of matrix A12, the second to the last column of Ab12.
q_of_uH = [-100] + [0 for nuH in range(N-2)] # the nodal demand q for unknown heads. For Node 1 is -100, others 0.
q_of_uH = scipy.sparse.csr_matrix(q_of_uH).T.tocsr()

nuQ = p                            # the number unknown flows. All flows are assumed unknown

A21 = A12.transpose()              #A21 is the transpose of A12
AH = A10*scipy.sparse.csr_matrix(H0) #AH is the A10H0 matrix composition of A10 and H0


def myF(z,a,r,nuQ,nuH,A12,A21,AH,q_of_uH):                        #function that will be given to the solver with output an array of equations

    v_Q = []                       # list of unknown pipe discharges, as the Q in Eq 16 in the Todini paper
    v_H = []                       # list of unknown nodal heads, as the H in Eq 16 in the Todini paper
    v_A11 = []                     # diagonal elements in the matrix A11    
    
    for j in range(0, nuQ):
        v_A11.append(r*pow(z[j],a-1))   # calculate the diagonal elements in A11
        v_Q.append(z[j])                   # update the pipe discharges based on the inital guess/iteration results
    
    for i in range(0, nuH):
        v_H.append(z[i+nuQ])              # update the nodal heads based on the inital guess/iteration results

    A11 = scipy.sparse.diags(v_A11)           # create the A11 diagonal matrix based on diagonal elements
    A1112 = scipy.sparse.hstack([A11,A12])        # combine A11 with A12, as in Eq 16 in the Todini paper

    dim_AZ = (nuH,nuH)             # dimensions of the lower right all-zero matrix in Eq 16
    AZ = scipy.sparse.csr_matrix(dim_AZ)          # create the lower right all-zero matrix
    A21Z = scipy.sparse.hstack([A21,AZ])

    K = scipy.sparse.vstack([A1112, A21Z]).tocsr() # the LHS matrix part in Eq 16

    AHq = scipy.sparse.vstack([AH,q_of_uH]).tocsr()    # the RHS of Eq 16
    
    QH = np.append(v_Q,v_H)        # update the unknowns in Eq 16 based on initial guess or iteration results, consisting of unknown flows and unknown heads
    QH = scipy.sparse.csr_matrix(QH).T.tocsr()
    
    KQH = K*QH          # the LHS of Eq 16
    
    F = KQH+AHq          # LHS - RHS, should equal to zeros

    F_lst = F.toarray().reshape((-1,))

    # print(KQH)
    # print(AHq)
    # sys.exit(0)
    
    return F_lst # Change to generalize


zGuess = ([100]*(p+N-1))          # initial guess of the unknowns, with only the head of one node is known
z = fsolve(myF, zGuess, args=(a,r,nuQ,nuH,A12,A21,AH,q_of_uH))

print(z) 
