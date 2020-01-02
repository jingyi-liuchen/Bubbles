#ifndef FIELD_H
#define FIELD_H

#include<cstdlib>
#include<cstdio>

#include "frame.h"

class Field
{
    public:
    int nmesh, fieldmesh[3];
    double fieldlen[3], meshsize[3];
    int** fieldnum;
    double  **fieldv, *fieldT;
 
    int nmesh_all; 
    int** fieldnum_all;                        
    double **fieldv_all, *fieldT_all;

    class Frame* frame;

    Field(class Frame* iframe);
    ~Field();


    void cal_field(FILE* outfield);   //driver function
    void cal_field_local();
    void comm_field();                //communicate 
    void cleanup_field();
    void output_field(FILE* outfield);

    private:
    Field(const Field&);
};

#endif
