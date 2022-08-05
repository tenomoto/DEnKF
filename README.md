# DEnKF

This repository contains the linear advection (LA) model and the quasi-geostropic (QG) model used with the deterministic ensemble Kalman filter (DEnKF) in Sakov and Oke (2008).
The code is used in the 26-th data assimilation (DA) summer school in 2022 held in Mutsu, Aomori, Japan and online.

## Models and DA algorithm

The models and assimilation code are mainly written in Python
with QG model also in Fortran. 

- Linear advection model (Evensen 2004)
- Quasi-Geostrophic model

## Files

- denkf.py: DEnKF 
- la.py: LA model
- plotstate.py: plots states of truth, initial guess, analysis, and observations.
- ploterror.py: plots time evolution of L2 and ensemble spread
- qg.py: QG model
- fd.py: five-point finite difference for Laplacian and Arakawa Jacobian
- ode.py: fourth-order Runge-Kutta
- mg.py: multigrid solver
- plotqg.py: plots vorticity and stream function
- fortran/qg_module.f90: components of QG models
- fortran/fd_module.f90: five-point finite difference for Laplacian and Arakawa Jacobian
- fortran/ode_module.f90: fourth-order Runge-Kutta
- fortran/mg_module.f90: multigrid solver
- fortran/run_qg.f90: driver for QG
- fortran/test_mg.f90: test program for multigrid

## References

-
