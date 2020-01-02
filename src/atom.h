#ifndef ATOM_H
#define ATOM_H

#include "frame.h"

class Atom
{
    public:
    int natoms;               //total number of atoms in this system
    int nlocal;               //number of local atoms
    int nremote;              //number of remote atoms
    int nremote_max;          //max number of remote atoms
    int nlocal_max;           //max number of local atoms

    class Frame* frame;
    int *id;                  //id in original lammpstrj
    int *type;                //atom type
    double **x, **v, **f;
    double **xre;             //remote atoms
 
    Atom(class Frame* iframe): frame{iframe},id{NULL},type{NULL}, x{NULL}, v{NULL},f{NULL},xre{NULL} {}
    ~Atom();
    
    void init();
    void add_atom_local(int aid, int atype, double* coord, double* vel, double* force);
    void init_remote(int* ncellx, int* nlayer, int* procgrid);
    void add_atom_remote(double* coord);
    void grow_nmax_chunk_local();
    void cleanup_atom();
    void print_atom_local();

    private:
    Atom(const Atom& );
};

#endif
