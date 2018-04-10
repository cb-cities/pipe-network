#fixed topology, variable unknown 
#TODO: take out the known heads, pass parameters to the function in stead of refering to the global value
#TODO: generalize the numbers based on the numbers of unknows

import numpy as np
from scipy.optimize import fsolve
import sys

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
print 'uH', uH
print 'A12', A12
print 'A10', A10
print 'H0', H0

print 'q_of_uH', q_of_uH
    
uQ = []                            #numbers of the unknown flows
nuQ = input('Insert the total number of unknown flows and press enter: ')   #total Number of Unknown Flows

for i in range(0, nuQ):
    c = input('Insert the pipe number which flow is unknown and press enter: ')
    uQ.append(c)
print 'uQ', uQ

A21 = A12.transpose()              #A21 is the transpose of A12
print 'A21', A21
AH = np.dot(A10,H0)                #AH is the A10H0 matrix composition of A10 and H0
print 'AH', AH

def myF(z):                        #function that will be given to the solver with output an array of equations
    
    print 'z', z
    v_Q = Q                        #
    v_H = []                       #
    v_A11 = []                     #array with the values of the diagonal in the matrix A11
    #w = 0                          #total number of unknown     
    
    for i in range(0, len(Q)):
        v_A11.append(r[i]*pow(Q[i],a-1)) #Bing
    for j in uQ:
        v_A11[j] = r[j]*pow(z[j],a-1)
        v_Q[j] = z[j]
        #w+=1
        # for j in uQ:
        #     if i == j:
        #         v_A11.append(r[i]*pow(z[w],a-1))
        #         v_Q[j] = z[w]
        #         w +=1
        #         print w
        #         break
        # else:
        #     v_A11.append(r[i]*pow(Q[i],a-1))
    
    for i in range(0, len(uH)):
        v_H.append(z[i+len(Q)]) ## Change to generalize, after taking out knows, change to uQ
        #w +=1
        #print '2nd', w
    #print 'v_A11', v_A11
    #print 'v_Q', v_Q
    #print 'v_H', v_H

    A11 = np.diag(v_A11)
    #print 'A11', A11
    A1112 = np.append(A11,A12,axis=1)
    #print 'A1112', A1112
    #A1112 = np.around(A1112,decimals=2)
    dim_AZ = (nuH,nuH)
    AZ = np.zeros(dim_AZ)
    #print 'AZ', AZ
    A21Z = np.append(A21,AZ,axis=1)
    #print 'A21Z', A21Z
    K = np.concatenate((A1112, A21Z), axis=0)
    #print 'K', K
    #sys.exit(0)

    AHq = np.append(AH,q_of_uH)
    t_AHq = AHq[np.newaxis, :].T
    #print 'AHq', AHq
    #print 't_AHq', t_AHq
    
    QH = np.append(v_Q,v_H)
    t_QH = QH[np.newaxis, :].T
    #print 'QH', QH
    #print 't_QH', t_QH
    
    KQH = np.matmul(K,t_QH)
    #print 'KQH', KQH
    
    F = np.add(KQH,t_AHq)
    t_F = F.T
    #print 'F', F
    #sys.exit(0)
    
    return t_F.reshape(len(z)) # Change to generalize

#zGuess = ([79,19,11,39,19,11,10,97])
zGuess = ([100]*(len(Q)+len(uH))) # Change to generalize
#print zGuess.shape
print zGuess
z = fsolve(myF, zGuess)

print z 
