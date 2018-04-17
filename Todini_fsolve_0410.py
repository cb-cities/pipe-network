#!Python2
### Based on the Global Algorithm Representation in Todini and Rossman (2013) (Eq 16 in the paper).
### According to the paper, for a water pipeline network with N nodes and p edges, no flow is observed. All pipe flows are unknown (the unknown flow's pipe index number stored in list varaible uQ). At least one nodal head is known. The index of the unknown heads are stored in the list variable uH. The number of unknown heads (nuH) satisfies 0 < nuH < N.
### The number of equations (mass conservation and energy conservation) equal to the number of unknowns (nuQ + nuH). Thus the system of equations can be solved using scipy.optimize.fsolve.

import numpy as np
from scipy.optimize import fsolve
import sys

a = 1.852       #coeff Hazen-Williams
r = [0.000298867,0.007789731,0.014060475,0.002157818,   #coefficient that depends on the dimensions used as well
     0.011684596,0.02812095,0.014060475]                #as on the pipe diameter, length, and possibly the flow rate
Q = [80, 20, 10, 40, 20, 10, 10]   #true flows per each pipe
q = [100, -10, -20, -30, -40]      #demand per each node
H = [100, 99, 98, 97, 96]          #hydraulic heads per each node
                        
p = 7                              # total number of pipes
N = 5                              # total number of nodes
Ab12 = np.array([                  #fixed topology, p x N matrix
        [-1, 1, 0, 0, 0],
        [-1, 0, 1, 0, 0],
        [ 0,-1, 1, 0, 0],
        [ 0,-1, 0, 1, 0],
        [ 0,-1, 0, 0, 1],
        [ 0, 0,-1, 0, 1],
        [ 0, 0, 0,-1, 1]])
                                     
A10 = Ab12                         #creation of matrix A10 from Ab12
H0 = H                             #creation of the vector H0 from H

nuH = input('Insert the total number of unknown heads (between 1 and {}) and press enter: '.format(N-1))  #total Number of Unknown Head
A12 = np.zeros((p, nuH))

uH = []                            # node numbers of the unknown head
q_of_uH = []                       # list of known demand (when nodal head is unknown, then the demand must be known)

for i in range(0, nuH): 
    c = input('Insert the node number of the unknown heads and press enter, index starting from 0: ')
    uH.append(c)
                                                
    A12[:, i] = Ab12[:,c]          #substitute sequentially the columns in the matrix A12 with the corresponding column of the unknown head 
    A10 = np.delete(A10, np.s_[c-i],1)    #A10 created by deleting sequentially the corresponding column of the unknown head
    H0.pop(c-i)                    #remove at data at position c-i
    
    q_of_uH.append(q[c])           # update list of known demand

### Bing modified this section to make all flows unknown, original code is commented out
# uQ = []                            #numbers of the unknown flows
# nuQ = input('Insert the total number of unknown flows and press enter: ')   #total Number of Unknown Flows

# for i in range(0, nuQ):
#     c = input('Insert the pipe number which flow is unknown and press enter: ')
#     uQ.append(c)
# print 'uQ', uQ
nuQ = p                            # the number unknown flows, basically equals to the number of pipes
uQ = [i for i in range(nuQ)]       # the list of pipe index for unknown flows, basically all pipes

A21 = A12.transpose()              #A21 is the transpose of A12
AH = np.dot(A10,H0)                #AH is the A10H0 matrix composition of A10 and H0

def myF(z,a,r,nuQ,nuH,A12,A21,AH,q_of_uH):                        #function that will be given to the solver with output an array of equations
    
    #print 'z', z
    v_Q = []                       # list of unknown pipe discharges, as the Q in Eq 16 in the Todini paper
    v_H = []                       # list of unknown nodal heads, as the H in Eq 16 in the Todini paper
    v_A11 = []                     # diagonal elements in the matrix A11    
    
    for j in range(0, nuQ):
        v_A11.append(r[j]*pow(z[j],a-1))   # calculate the diagonal elements in A11
        v_Q.append(z[j])                   # update the pipe discharges based on the inital guess/iteration results
    
    for i in range(0, nuH):
        v_H.append(z[i+nuQ])              # update the nodal heads based on the inital guess/iteration results

    A11 = np.diag(v_A11)           # create the A11 diagonal matrix based on diagonal elements
    A1112 = np.append(A11,A12,axis=1)        # combine A11 with A12, as in Eq 16 in the Todini paper

    dim_AZ = (nuH,nuH)             # dimensions of the lower right all-zero matrix in Eq 16
    AZ = np.zeros(dim_AZ)          # create the lower right all-zero matrix
    A21Z = np.append(A21,AZ,axis=1)

    K = np.concatenate((A1112, A21Z), axis=0) # the LHS matrix part in Eq 16

    AHq = np.append(AH,q_of_uH)    # the RHS of Eq 16
    t_AHq = AHq[np.newaxis, :].T
    
    QH = np.append(v_Q,v_H)        # update the unknowns in Eq 16 based on initial guess or iteration results, consisting of unknown flows and unknown heads
    t_QH = QH[np.newaxis, :].T
    
    KQH = np.matmul(K,t_QH)        # the LHS of Eq 16
    
    F = np.add(KQH,t_AHq)          # LHS - RHS, should equal to zeros
    t_F = F.T
    
    return t_F.reshape(len(z)) # Change to generalize


zGuess = ([100]*(nuQ+nuH))          # initial guess of the unknowns
#print 'Initial Guess {}'.format(zGuess)
z = fsolve(myF, zGuess, args=(a,r,nuQ,nuH,A12,A21,AH,q_of_uH))

print z 
