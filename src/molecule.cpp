#include "molecule.h"
#include "function.h"

Molecule::~Molecule()
{
    if(beadmass) free_1D<double>(beadmass);
    if(atomtype) free_1D<int>(atomtype);
    if(chemid)   free_1D<int>(chemid);
    if(atomname) free_2D<char>(atomname);
}

void Molecule::init()
{
    allocate_1D<double>(beadmass, nunit, "molecule->beadmass");
    allocate_1D<int>(atomtype, nunit, "molecule->atomtype");
    allocate_1D<int>(chemid, nunit, "molecule->chemid");
    allocate_2D<char>(atomname,nunit,8,"molecule->atomname");
}
