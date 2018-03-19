#fixed topology, variables unknown 

import numpy as np
from scipy.optimize import fsolve

a = 1.852       #coeff Hazen-Williams
r = [0.000298867,0.007789731,0.014060475,0.002157818,   #coefficient that depends on the dimensions used as well
     0.011684596,0.02812095,0.014060475]                #as on the pipe diameter, length, and possibly the flow rate
Q = [80, 20, 10, 40, 20, 10, 10]   #true flows per each pipe
q = [100, -10, -20, -30, -40]      #demand per each node
H = [100, 99, 98, 97, 96]          #hydraulic heads per each node

Ab12 = np.array([
        [-1, 1, 0, 0, 0],
        [-1, 0, 1, 0, 0],
        [ 0,-1, 1, 0, 0],
        [ 0,-1, 0, 1, 0],
        [ 0,-1, 0, 0, 1],
        [ 0, 0,-1, 0, 1],
        [ 0, 0, 0,-1, 1]])

A10 = Ab12
H0 = H

nuH = input('Insert the total number of unknown heads and press enter: ')
A12 = np.zeros((len(Ab12), nuH))

uH = []

for i in range(0, nuH): 
    c = input('Insert the number of the unknown head and press enter: ')
    uH.append(c)
                                                
    A12[:, i] = Ab12[:,c]
    A10 = np.delete(A10, np.s_[c-i],1)
    H0.pop(c-i)
    #print Ab12
    #print A12
    #print A10
    #print uH
    #print len(uH)
    print H0
    
uQ = []
nuQ = input('Insert the total number of unknown flows and press enter: ')

for i in range(0, nuQ):
    c = input('Insert the pipe number which flow is unknown and press enter: ')
    uQ.append(c)
    #print uQ

A21 = A12.transpose()
AH = np.dot(A10,H0)
print AH
