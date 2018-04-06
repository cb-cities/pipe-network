#fixed topology, variable unknown 

import numpy as np
from scipy.optimize import fsolve

a = 1.852       #coeff Hazen-Williams
r = [0.000298867,0.007789731,0.014060475,0.002157818,   #coefficient that depends on the dimensions used as well
     0.011684596,0.02812095,0.014060475]                #as on the pipe diameter, length, and possibly the flow rate
Q = [80, 20, 10, 40, 20, 10, 10]   #true flows per each pipe
q = [100, -10, -20, -30, -40]      #demand per each node
H = [100, 99, 98, 97, 96]          #hydraulic heads per each node
                        
Ab12 = np.array([                  #fixed topology
        [-1, 1, 0, 0, 0],
        [-1, 0, 1, 0, 0],
        [ 0,-1, 1, 0, 0],
        [ 0,-1, 0, 1, 0],
        [ 0,-1, 0, 0, 1],
        [ 0, 0,-1, 0, 1],
        [ 0, 0, 0,-1, 1]])
                                     
A10 = Ab12                         #creation of matrix A10 from Ab12
H0 = H                             #creation of the vector H0 from H

nuH = input('Insert the total number of unknown heads and press enter: ')  #total Number of Unknown Head
A12 = np.zeros((len(Ab12), nuH))

uH = []                            #numbers of the unknown head
q_of_uH = []                       #q related to the unknown heads

for i in range(0, nuH): 
    c = input('Insert the number of the unknown head and press enter: ')
    uH.append(c)
                                                
    A12[:, i] = Ab12[:,c]          #substitute sequentially the columns in the matrix A12 with the corresponding column of the unknown head 
    A10 = np.delete(A10, np.s_[c-i],1)    #A10 created by deleting sequentially the corresponding column of the unknown head
    H0.pop(c-i)                    #remove at data at position c-i
    
    q_of_uH.append(q[c])
    
uQ = []                            #numbers of the unknown flows
nuQ = input('Insert the total number of unknown flows and press enter: ')   #total Number of Unknown Flows

for i in range(0, nuQ):
    c = input('Insert the pipe number which flow is unknown and press enter: ')
    uQ.append(c)

A21 = A12.transpose()              #A21 is the transpose of A12
AH = np.dot(A10,H0)                #AH is the A10H0 matrix composition of A10 and H0


def myF(z):                        #function that will be given to the solver with output an array of equations
    
    v_Q = Q                        #
    v_H = []                       #
    v_A11 = []                     #array with the values of the diagonal in the matrix A11
    w = 0                          #total number of unknown     
    
    for i in range(0, len(Q)):
        for j in uQ:
            if i == j:
                v_A11.append(r[i]*pow(z[w],a-1))
                v_Q[j] = z[w]
                w =+1
                break
        else:
            v_A11.append(r[i]*pow(Q[i],a-1))
    
    for i in range(0, len(uH)):
        v_H.append(z[w])
        w =+1

    A11 = np.diag(v_A11)
    A1112 = np.append(A11,A12,axis=1)
    #A1112 = np.around(A1112,decimals=2)
    dim_AZ = (nuH,nuH)
    AZ = np.zeros(dim_AZ)
    A21Z = np.append(A21,AZ,axis=1)
    K = np.concatenate((A1112, A21Z), axis=0)

    AHq = np.append(AH,q_of_uH)
    t_AHq = AHq[np.newaxis, :].T
    
    QH = np.append(v_Q,v_H)
    t_QH = QH[np.newaxis, :].T
    
    KQH = np.matmul(K,t_QH)
    
    F = np.subtract(KQH,t_AHq)
    t_F = F.T
    
    return t_Fa

zGuess = ([100,100,100,100,100,100,100,100,100])
#print zGuess.shape
print zGuess
z = fsolve(myF,zGuess)

print z 
