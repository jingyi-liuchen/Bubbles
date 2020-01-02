#ifndef FUNCTION_H
#define FUNCTION_H

#include <cstdio>
#include <cstdlib>
#include<cinttypes>
#include<cstdint>

void check_bound(double* coor,const double* blo, const double* bhi);
void  mimage(double* delx, const double* boxlen);
double cal_distsq(const double* coor1,const double* coor2, const double* boxlen);
int mapto1D(const int* id_3D, const int* mesh);
void mapto3D(int id_1D, int* id_3D, const int* mesh);
void gridtocms(const int* id_3D, const double* meshsize, const double* sublo, double* cms);
void id_pbc(int& xid, int nx);

void* smalloc(int64_t nbytes, const char* name);
void* srealloc(void* ptr, int64_t nbytes, const char* name);

template <typename T>
T* allocate_1D(T*& array, int size1, const char* name)
{
    int64_t nbytes = sizeof(T)*size1;
    array = (T*) smalloc(nbytes,name);
    return array;
}

template <typename T>
void free_1D(T*& array)
{
    free(array);
    array = NULL;
}

template <typename T>
T* grow_1D(T*& array, int size1, const char* name)
{
    if(array==NULL)
        return allocate_1D(array, size1, name);

    int64_t nbytes = sizeof(T)*size1;
    array = (T*) srealloc(array,nbytes,name);
    return array;
}

template <typename T>
T** allocate_2D(T**& array, int size1,int size2,const char* name)
{
    int64_t nbytes = sizeof(T)*size1*size2;
    T* data  = (T*) smalloc(nbytes, name);
    nbytes = sizeof(T*)*size1;
    array = (T**) smalloc(nbytes, name);

    int64_t n=0;
    for(int i=0;i<size1;i++)
    {
        array[i] = &data[n];
        n += size2;
    }
    return array;
}

template <typename T>
void free_2D(T**& array)
{
    free(array[0]);
    free(array);
    array = NULL;
}

template <typename T>
T** grow_2D(T**& array, int size1, int size2, const char* name)
{
    if (array == NULL) return allocate_2D(array,size1,size2, name);

    int64_t nbytes =  sizeof(T) * size1* size2;
    T *data = (T*) srealloc(array[0],nbytes,name);
    nbytes = sizeof(T*) * size1;
    array = (T**) srealloc(array,nbytes,name);

    int64_t n=0;
    for(int i=0;i<size1;i++)
    {
        array[i] = &data[n];
        n += size2;
    }
    return array;
}

#endif
