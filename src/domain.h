#ifndef DOMAIN_H
#define DOMAIN_H

#include "frame.h"

class Domain
{
    public:
    double boxlo[3], boxhi[3], sublo[3], subhi[3], boxlen[3], sublen[3];  
    class Frame* frame;

    Domain(class Frame* iframe); 
 
    private:
    Domain(const Domain &); 
};



#endif
