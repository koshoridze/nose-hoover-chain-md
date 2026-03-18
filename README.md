# Nose–Hoover Chain MD

A lightweight C++ molecular dynamics (MD) implementation featuring a Nose–Hoover chain thermostat and barostat.

## Overview

This project implements a simple molecular dynamics engine with:

* Nose–Hoover chain (NHC) thermostat
* Nose–Hoover barostat (NPT ensemble)
* Embedded Atom Method (EAM) potentials
* XYZ trajectory export

---

## Dependencies

* C++20 compatible compiler
* CMake ≥ 3.16
* Eigen3
* GSL

---

## Build

```bash
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

---

## Usage

After building, run:

```bash
./nosehoover-chain
```

(Adjust input parameters in the source if needed.)

---


## Notes

* This is a minimal implementation and not intended for production simulations
* Numerical stability and parameter tuning are left to the user
* Results should be validated before scientific use

## References

The integration scheme is based on:

1. Martyna, G. J., Tuckerman, M. E., Tobias, D. J., & Klein, M. L.  
   *Explicit reversible integrators for extended systems dynamics.*  
   Molecular Physics, 87(5), 1117–1157 (1996).  
   https://doi.org/10.1080/00268979600100761
