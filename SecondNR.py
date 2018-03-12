
# coding: utf-8

# In[8]:


from numpy import *
from scipy.optimize import fsolve
a = 1.852
r = [0.000298867,0.007789731,0.014060475,0.002157818,0.011684596,0.02812095,0.014060475]
AH = [1,2,1,99,3,2,-96]
A21 = [0,0,0,1,0,0,-1]


def myFunction(uQ):
    
    Q1 = uQ[0]
    Q2 = uQ[1]
    Q3 = uQ[2]
    Q4 = uQ[3]
    Q5 = uQ[4]
    Q6 = uQ[5]
 
    H3 = uQ[6]
 
    F=empty((7))
 
    F[6] = Q4-40
 
    for i in range(0,6):
        F[i] = r[i]*pow(uQ[i],a)-AH[i]
        if A21[i] != 0:
            F[i] = F[i]+A21[i]*H3 

    return F

zGuess = array([100,100,100,100,100,100,100])
z = fsolve(myFunction,zGuess)
print(z)

