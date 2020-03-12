#ifndef FRAME_H
#define FRAME_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <mpi.h>

#include "input.h"
#include "atom.h"
#include "domain.h"
#include "bubble.h"
#include "shell.h"
#include "field.h"

class Frame
{
    public:
    int nproc;                                        //number of procs             
    int myid;                                         //my rank in MPI_COMM_WORLD 
    MPI_Comm cart;                                    //communicator in cartesian topology
    int procgrid[3];                                  //number of procs in each dimension
    int myid_3D, mycor[3];                            //my rank in cartesian

    FILE* intraj;                                     //trajetory file

    int index;                                        //current timestep
    char inpara_file[32];                             //input parameter file

    class Input*  input;                              //input file parameters
    class Domain* domain;                             //domain and box size information
    class Atom* atom;                                 //atom-related information

    class Bubble* bubble;                             //calculate volume, cms, and kappa for bubble(s)
    class Shell* shell;                               //calculate shell number, temperature, radial velicity 
    class Field* field;                               //calculate field data of number, temperature, velocity
   
    Frame(const char* inpara_file);
    ~Frame();

    //allocate memory for pointers  and set up cartesian communicator
    void init();

    //set up processor 3D grid
    void set_procgrid();
    int factor(int nproc, int**factors);
    int best_factors(int npossible, int **factors, int *best);

    //read trajectory file and set local atom information
    int read_head();
    int read_atom();
    void set_atom_info(int nline, char* buf);

    private:
    Frame(const Frame &);                            //prohibit use of copy constructor
};

#endif
