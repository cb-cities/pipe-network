from scipy.optimize import fsolve
import math

def equations(p):
    x, y = p
    return (x+y-4, x*y - 3)

x, y =  fsolve(equations, (1, 1))

print equations((x, y)), x, y