#ifndef SHELL_H
#define SHELL_H

#include <cstdlib>
#include <cstdio>
#include<map>

class Shell
{
    public:
    int nshell, del, n_inner_shell, natom_type;
    std::map<int, int> atom_type_map;
    std::map<int, double> atom_mass_map;
    double cms[3];
    double* shellr, *shellT, *shellvr;
    int** shellnum;
    class Frame* frame;

    Shell(class Frame* iframe);
    ~Shell();

    void cal_shell(FILE* outshell);        //driver function
    void cal_shellr();                         
    int cal_shellid(double dist);
    void cal_shellprop();
    void cal_shellr_deln();
    void cleanup_shell();
    void output_shell(FILE* outshell);

    private:
    Shell(const Shell&);
};

#endif
