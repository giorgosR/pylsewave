# PylseWave

![pylsewavelogo](./report/DocOnce/figures/pylsewave_logo.png =50x50)

### A C++/Python package to solve pulse wave dynamics in arterial networks

The following diagram depicts the rationale along with the structure of the package

![pylsewave toolkit](./JupyterNbs/images/pylsewave.png)

A c/python library to solve 1D pulse wave propagation in blood vessels or any hyperbolic PDE system in the form of

![Hyperbolic system](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20%5Cbf%7BU%7D%7D%7B%5Cpartial%20t%7D%20&plus;%20%5Cfrac%7B%5Cpartial%20%5Cbf%7BF%7D%7D%7B%5Cpartial%20x%7D%20%3D%20%5Cbf%7BS%7D)

### Purpose
This library is designed to solve efficiently pulse wave dynamics in human arteries. It is written in python, cython and there are some c++ classes.

### Solver Example
```python
# ---- SOLUTION WITH CYTHON CLASSES ----- #
from pulsewavepy.cynum import cPDEsWat, BCsADAN56, cMacCormackSolver

myPDEs = cPDEsWat(vesssel_network)
myBCs = BCsADAN56(myPDEs, q_inlet_bc.eval_spline)

U0_vessel = np.array([0],dtype=np.int)
UL_vessel = np.array(terminal_vessels.keys())
UBif_vessel = np.array(bif_vessels)
UConj_vessel = np.array(conj_points)

mySolver = cMacCormackSolver(myBCs)
mySolver.set_T(dt=dt, T=T, no_cycles=tc)
mySolver.set_BC(U0_vessel, UL_vessel, UBif_vessel, UConj_vessel)
# ------- SOLVE AND TIME --------------- #
%time mySolver.solve(casename, myCallback, cfl_n=0.6, stability_safety_factor=1.0)
myCallback.close_file(casename)
```