# Project Title
This program reads in lammps trajectory file and (optionally) output the following:   
(1) number, size, and shape anisotropy factor for bubbles (void clusters);  
(2) shell density, temperature, and radial velocity profile;  
(3) 3D field data of density, temperature and velocity.  

Details regarding analysis methods can be found in   
Chen, Jingyi L., et al. "Molecular Simulations Probing the Thermophysical Properties of Homogeneously Stretched and Bubbly Water Systems." Journal of Chemical & Engineering Data (2019).

### Prerequisites

MPI and Eigen library is required for compiling.

### Installing

Modify the first two lines in Makefile to specifiy c++ compilers and path of Eigen library.  
Type "make" for compiling.  
 
## Running the tests

Type "make test" for automated test.
