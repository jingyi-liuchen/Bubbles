#ifndef BUBBLE_H
#define BUBBLE_H

#include "frame.h"

class Bubble
{
    public:
    int ncellx[3], ncell;                             //number of link cell and its 3D topology
    int* icell, *nicell;                             //cell index for moleules, number of molecules in a cell
    int** iucell;                                     //iucell(x,y) mol id for yth molecule in cell x 
    int nlayer[3];                                    //number of layers for surface cells

    int* allcellty;                                   //cell type (only allocated for rank 0)
    int ncluster,realncut;                            //number of bubbles and actually outputted bubbles
    std::vector<int> vcluster;                        //volume of bubbles
    std::vector<std::vector<int>> clustermem;         //member of cluster
    double* kappa;                                    //shape anisotropy factor for all bubbles
    double** ccms;                                    //center of mass for all bubbles

    double maxcms[3];                                    //center of mass for biggest bubble, used for later shell calculation 
    int maxvol;                                       //volume for largest bubble  

    class Frame* frame;

    Bubble(class Frame* iframe);
    ~Bubble();

    void cal_bubble(FILE* outbubble);                                    //driver function
    void BFS();                                                          //breath first search (driver)
    int BFS_visit(int id_1D, const int* mesh, int* visited);
    int* sort_vol();                                                     //sort bubble based on volume
    void cal_cms(const int* sortid);                                     //calculate bubble cms 
    void cal_kappa(const int* sortid);                                   //calculate bubble shape anisotropy factor
    void output_vol(FILE* outbubble);                                    //output bubble info  
    void cleanup_bubble();                                               //clean up memory for bubble calculation

    void set_linkcell();                                                 //set up link cell for calculating nearest neighbor
    void comm_surfz();
    void comm_surfy();
    void comm_surfx();

    int cal_molty(int imol, int& id_1D);                               //decide if a molecule is vapor or liquid
    void cal_cellty();                                                 //decide if a mesh is vapor or liquid
    void comm_cellty(const int* cellty);
    int comm_surfz_molty(const int* molty, int* remolty, int offset);
    int comm_surfy_molty(const int* molty, int* remolty, int offset);
    int comm_surfx_molty(const int* molty, int* remolty, int offset);

    private:
    Bubble(const Bubble&);
};




#endif
