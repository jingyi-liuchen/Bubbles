#ifndef INPUT_H
#define INPUT_H

#include<cstdio>
#include<cstdlib>
#include "molecule.h"
#include "frame.h"

class File_para
{
    public:
    char traj_file[32];
    int maxframe;

    File_para():maxframe{1} {}
};

class Bubble_para
{
    public:
    int lbubble[3];
    int bubble_mesh[3],ncut,nliqcut;
    double rcut;

    Bubble_para():lbubble{1,0,0},bubble_mesh{20,20,20},ncut{1},
                  nliqcut{2},rcut{3.3} {}
};

class Shell_para
{
    public:    
    int lshell[3],del;
    int nshell;
    double delr,delv;
    int deln;

    Shell_para():lshell{1,0,0},del{0},nshell{20},delr{0.2},delv{100.0},deln{1000} {}
};

class Field_para
{
    public:
    int lfield[3];
    int field_mesh[3];

    Field_para():lfield{1,0,0},field_mesh{4,4,4} {}
};

class Sys_para
{
    public:
    int init_flag;
    int nmolty, nmolty_read, natoms;
    double boxl[3];
    Molecule*  all_molecule;

    Sys_para():init_flag{0},nmolty{1},nmolty_read{0}, natoms{0},
               boxl{40.0,40.0,40.0}, all_molecule{NULL} {}
    ~Sys_para(); 
    void init();

    private:
    Sys_para(const Sys_para&);
};

class Input
{
    public:
    char          inpara_file[32];
    FILE*         inpara;
    class Frame*        frame;
    File_para*    file_pars;
    Bubble_para*  bubble_pars; 
    Shell_para*   shell_pars;
    Field_para*   field_pars;
    Sys_para*     sys_pars;

    Input(Frame* iframe):inpara{NULL},frame{iframe},file_pars{NULL},bubble_pars{NULL},
                        shell_pars{NULL},field_pars{NULL},sys_pars{NULL} {}
    ~Input();
   
    void init();
    int read_section(const char* sec_name);
    void set_sec_para(const char* sec_name,  char* lines);
    void set_line_para(const char* sec_name, char tokens[16][32], int ntokens);
    void print_message(int read_result, const char* sec_name);
    void strip_newline(const char* oldline, char* newline);
    void print_tokens(char tokens[16][32], int ntokens);

    private:
    Input(const Input&);
};

#endif
