#include "atom.h"
#include "function.h"

#define LB_FACTOR 1.1
#define DELTA 16384
#define MAXINCELL 5

Atom::~Atom()
{
    if(id)    free_1D<int>(id);
    if(type)  free_1D<int>(type);
    if(x)     free_2D<double>(x);
    if(v)     free_2D<double>(v);
    if(f)     free_2D<double>(f); 
    if(xre)   free_2D<double>(xre);
}

void Atom::init()
{
    natoms = frame->input->sys_pars->natoms;
    int nproc = frame->nproc;

    nlocal = 0;
    nremote = 0;

    if(nproc==1)
    {
        nlocal_max = natoms;
    }
    else
    {
        nlocal_max = natoms/nproc * LB_FACTOR;
    }

    allocate_1D<int>(id,nlocal_max,"atom->id");
    allocate_1D<int>(type,nlocal_max,"atom->type");
    allocate_2D<double>(x, nlocal_max, 3, "atom->x");
    allocate_2D<double>(v, nlocal_max, 3, "atom->v");
    allocate_2D<double>(f, nlocal_max, 3, "atom->f");
}

void Atom::add_atom_local(int aid, int atype, double* coord, double*vel, double* force)
{
    id[nlocal] = aid;
    type[nlocal] = atype;
    for(int i=0;i<3;i++)
    {
        x[nlocal][i] = coord[i];
        v[nlocal][i] = vel[i];
        f[nlocal][i] = force[i];
    }
    nlocal ++;
    if(nlocal==nlocal_max)
        grow_nmax_chunk_local();  
}

void Atom::init_remote(int* ncellx, int* nlayer, int* procgrid)
{
   nremote_max = 0;
   for(int i=0;i<3;i++)
   {
      nremote_max += ncellx[(i+1)%3]*ncellx[(i+2)%3]*nlayer[i]*MAXINCELL*2;  
   }
   allocate_2D<double>(xre,nremote_max,3, "atom->xre");  
}

void Atom::add_atom_remote(double* coord)
{
    for(int i=0;i<3;i++)
    {
        xre[nremote][i] = coord[i];
    }
    nremote++; 
}

void Atom::grow_nmax_chunk_local()
{
  nlocal_max = nlocal_max/DELTA * DELTA;
  nlocal_max += DELTA;
  x = grow_2D<double>(x,nlocal_max,3, "atom->x");
  v = grow_2D<double>(v,nlocal_max,3, "atom->v"); 
  f = grow_2D<double>(f,nlocal_max,3, "atom->f");
  type = grow_1D<int>(type,nlocal_max, "atom->type");
  id = grow_1D<int>(id,nlocal_max, "atom->id");
}

void Atom::cleanup_atom()
{
    if(id)    free_1D<int>(id);
    if(type)  free_1D<int>(type);
    if(x)     free_2D<double>(x);
    if(v)     free_2D<double>(v);
    if(f)     free_2D<double>(f);
    if(xre)   free_2D<double>(xre);
}

void Atom::print_atom_local()
{
    char const * FORMAT = "%10d %10d %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n";
    printf("I have %d atoms\n",nlocal);
    for(int i=0;i<nlocal;i++)
    {
        printf(FORMAT,id[i],type[i],x[i][0],x[i][1],x[i][2],v[i][0],v[i][1],v[i][2],f[i][0],f[i][1],f[i][2]);
    }

}
