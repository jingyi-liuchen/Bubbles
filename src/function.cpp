#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "function.h"

int mapto1D(const int* id_3D, const int* mesh)
{
    return (id_3D[0]*mesh[1]*mesh[2] + id_3D[1]*mesh[2] + id_3D[2]);
}

void mapto3D(int id_1D, int* id_3D, const int* mesh)
{
    int nx = mesh[0], ny = mesh[1], nz = mesh[2];
    id_3D[0] = id_1D/(ny*nz);
    id_3D[1] = id_1D/nz - id_3D[0]*ny;
    id_3D[2] = id_1D%nz;
}

void mimage(double* delx, const double* boxlen)
{
    for(int i=0;i<3;i++)
    {
        delx[i] -= round(delx[i]/boxlen[i])*boxlen[i];
    }
}

double cal_distsq(const double* coor1,const double* coor2, const double* boxlen)
{
    double delx[3];
    for(int i=0;i<3;i++)
    {
        delx[i] = coor1[i] - coor2[i];
    }
    mimage(delx,boxlen);
    double distsq = 0.;
    for(int i=0;i<3;i++)
    {
        distsq += delx[i]*delx[i];
    }
    return distsq;
}

void gridtocms(const int* id_3D, const double* meshsize, const double* sublo, double* cms)
{
    for(int i=0;i<3;i++)
    {
        cms[i] = sublo[i] + id_3D[i] * meshsize[i] + 0.5*meshsize[i];
    }
}

void id_pbc(int& xid, int nx)
{
    while(xid < 0) {xid +=  nx;}
    while(xid > nx-1) {xid -= nx;}
}

void check_bound(double* coor, const double* blo, const double* bhi)
{
    for(int i=0;i<3;i++)
    {
        double len = bhi[i] - blo[i];
        if(coor[i] < blo[i])
        {
            coor[i] += len;
        }
        else if (coor[i] > bhi[i])
        {
            coor[i] -= len;
        }
    }
}

void* smalloc(int64_t nbytes, const char* name)
{
    if(nbytes==0) return NULL;
    void *ptr = malloc(nbytes);
    if(ptr==NULL)
    {
        printf("Failed to allocate %" PRId64 " bytes for %s",nbytes,name);
        exit(1);
    }

    return ptr;
}

void* srealloc(void* ptr, int64_t nbytes, const char* name)
{
    if(nbytes==0)
    {
        free(ptr);
        return NULL;
    }

    ptr = realloc(ptr,nbytes);
    if(ptr==NULL)
    {
        printf("Failed to allocate %" PRId64 " bytes for %s",nbytes,name);
        exit(1);
    }
    return ptr;
}

