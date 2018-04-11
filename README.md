# pipe-network
======

Solve for flows and unknown nodal heas in a water distribution network, based on [Todini & Rossman, 2013](Todini_2013.pdf).   

The current python script is [Todini_fsolve_0410.py](Todini_fsolve_0410.py). It utilize `scipy.optimize.fsolve` to solve the system of equations. It works with the general case when no pipe flow is observed, while at least one nodal head has to be known.   

When some pipe flows are observed, these known flows create redundancy in the systems of equations. Thus an target function needs to be supplied and optimization solvers such as `scipy.optimize.minimize` have to be used instead.
