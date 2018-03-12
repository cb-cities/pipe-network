
# coding: utf-8

# In[1]:


from numpy import *
from scipy.optimize import fsolve

def myFunction(z):
    q1=z[0]
    q2=z[1]
    q3=z[2]
    q4=z[3]
    q5=z[4]
    q6=z[5]
    h3=z[6]
    
    F=empty((7))
    F[0] = 0.00029887*pow(q1,1.852)-1
    F[1] = 0.00778973*pow(q2,1.852)-2
    F[2] = 0.01406048*pow(q3,1.852)-1
    F[3] = 0.00215782*pow(q4,1.852)+h3-99
    F[4] = 0.0116846*pow(q5,1.852)-3
    F[5] = 0.02812095*pow(q6,1.852)-2
    F[6] = q4-40
    return F

zGuess = array([100,100,100,100,100,100,100])
z=fsolve(myFunction,zGuess)
print(z)

