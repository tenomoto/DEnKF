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

- Evensen, G., 2004: Sampling strategies and square root analysis schemes for the EnKF. Ocean Dynamics, 54, 539–560, [doi:10.1007/s10236-004-0099-2](https://doi.org/10.1007/s10236-004-0099-2).
- Sakov, P., and P. R. Oke, 2008: A deterministic formulation of the ensemble Kalman filter: an alternative to ensemble square root filters. Tellus A, 60, 361–371, [doi:10.1111/j.1600-0870.2007.00299.x](https://doi.org/10.1111/j.1600-0870.2007.00299.x).

