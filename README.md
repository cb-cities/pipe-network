# pipe-network

Solve for flows and unknown nodal heas in a water distribution network, based on the Global Algorithm (Eq 16) in [Todini & Rossman, 2013](https://ascelibrary.org/doi/abs/10.1061/%28ASCE%29HY.1943-7900.0000703).   

The current python script is [Todini_fsolve_0410.py](Todini_fsolve_0410.py). It utilizes `scipy.optimize.fsolve` to solve the hydraulic system of equations representing the conservation of mass and the conservation of energy. It works with the case when no pipe flow is observed, while at least one nodal head has to be known.   

When some pipe flows are observed, these known flows create redundancy in the hyraulic system of equations. `scipy.optimize.fsolve` cannot handle over-determined systems of equations properly, thus a target function needs to be supplied and optimization solvers such as `scipy.optimize.minimize` have to be used instead.
