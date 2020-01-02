#ifndef MOLECULE_H
#define MOLECULE_H

#include<cstring>

class Molecule
{
    public:
    int init_flag,nunit, ncount,fix_dof;
    double* beadmass;
    int* atomtype, *chemid;
    char** atomname;

    Molecule():init_flag{0}, nunit{1},ncount{2000},fix_dof{0},beadmass{NULL},
               atomtype{NULL},chemid{NULL},atomname{NULL} {}
    ~Molecule();
    void init();

    private:
    Molecule(const Molecule&);
};

#endif
