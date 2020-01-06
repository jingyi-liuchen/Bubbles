## General description

There are five sections in input.dat: FILE, BUBBLE, SHELL, FIELD and SYSTEM. Comment line starts with #.    

### FILE

traj\_file: lammps trajectory file name  
maxfrmae: maximum number of frames to analyze  

### BUBBLE

lbubble: whether to calculate size, center of mass, shape anisotropy factor (0:no 1:yes)  
bubble\_mesh: number of meshes used in x, y, z direction  
ncut: only output information for ncut largest bubbles  
nliqcut: number of minimum nearest neighbor molecules to be considered as liquid molecule  
rcut: distance for defining nearest neighbor (ansgtrom)   

### SHELL

lshell: whether to calculate shell density, temperature, radial velocity profile  
del: 0: constant distance binning 1: constant volume binning 2: constant atom number binning  
nshell: number of shells  
delr: value of constant distance used when del = 0 (angstrom)  
delv: value of constant volume used when del = 1 (angstrom^3)  
deln: value of constant atom number used when del = 2  

### FIELD

lfield: whether to calculate field density, temperature, velocity  
field\_mesh: number of meshes used in x, y, z direction  

### SYSTEM

boxl: length of box in x, y, z directions (only used for deciding MPI processor grids, angstrom)  
nmolty: number of molecular types    

For each molecule type, specify one line starting with molty.   
nunit: number of atoms in a molecule  
ncount: number of molecules in the system  
beadmass: molar mass for each atom (g/mol) 
atomtype: type for each atom  
chemid: chemical id for each atom  
atomname: name for each atom  
