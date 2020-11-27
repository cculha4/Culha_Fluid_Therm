# Culha Fluid Thermodynamic Simulator

## Description

Repository for Culha et al., 2021 research on self-sustaining instability and its imprint on crystal zonations.
Project is current submitted to ESSOAr with the intention of submission to JGR. 


## Setup

Download the repository here.

## Run a simulation

In terminal you can run this script below:

```
gfortran -O4 2D-TwoPhase-GhostFluid-ST-FS-particle-gas-Cansu-vbc-enclave_1liq_wTracers_cntr_Unzen.f90 umf4_f77wrapper.o libumfpack.a libamd.a libsuitesparseconfig.a -lm -lrt
```

## Different simulations

You can compare the difference between the simulators, but attached you will find two fortran files. One has thermodynamics coupled:

```
2D-TwoPhase-GhostFluid-ST-FS-particle-gas-Cansu-vbc-enclave_1liq_wTracers_cntr_Unzen.f90
```

and the other does not: 

```
2D-TwoPhase-GhostFluid-ST-FS-particle-gas-Cansu-vbc-enclave_1liq.f90
```
