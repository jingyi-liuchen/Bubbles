# Project Title
This program Bubbles reads in lammps trajectory file and (optionally) output the following:   
(1) number, size, and shape anisotropy factor for bubbles (void clusters);  
(2) shell density, temperature, and radial velocity profile;  
(3) 3D field data of density, temperature and velocity.  

Details regarding analysis methods can be found in   
Chen, Jingyi L., et al. "Molecular Simulations Probing the Thermophysical Properties of Homogeneously Stretched and Bubbly Water Systems." Journal of Chemical & Engineering Data (2019).

## Prerequisites

MPI and Eigen library (http://eigen.tuxfamily.org) is required for compiling.

## Installing

Adapt Makefile to change c++ compilers (default:mpiicpc) and path of Eigen library (default:../)  

```
cd Bubbles
make   
```

## Running the tests

```
cd Bubbles
make test
```

## Example usage

```
mpirun -np 4 Bubbles/obj/bubble_ana input.dat
```  

## Description of input parameters
see Bubbles/tests/README.md
