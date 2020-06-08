#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <mpi.h>

#include "input.h"
#include "molecule.h"
#include "frame.h"

#define MAXLEN 1024

Sys_para::~Sys_para()
{
    if(all_molecule)  delete [] all_molecule;
}

void Sys_para::init()
{
    all_molecule = new Molecule[nmolty];
}

Input::~Input()
{
    if(file_pars)   delete file_pars;
    if(bubble_pars) delete bubble_pars;
    if(shell_pars)  delete shell_pars;
    if(field_pars)  delete field_pars;
    if(sys_pars)    delete sys_pars;
}

//optional sections:BUBBLE, SHELL, FIELD
void Input::init()
{
    strcpy(inpara_file,frame->inpara_file);
    int myid = frame->myid;

    if(myid==0)
    {
        inpara = fopen(inpara_file,"r");
        if(!inpara)
        {
            printf("Failed to open input parameter file!\n");
        }
    }

    int result;
    result = read_section("FILE");
    assert(result==0);

    result = read_section("BUBBLE");
    assert(result==0);

    result = read_section("SHELL");
    assert(result==0);

    result = read_section("FIELD");
    assert(result==0);

    result = read_section("SYSTEM");
    assert(result==0);
    assert(sys_pars->nmolty == sys_pars->nmolty_read);

    sys_pars->natoms = 0;
    Molecule* molecules  = sys_pars->all_molecule;
    for(int imolty=0; imolty<sys_pars->nmolty; imolty++)
    {
        sys_pars->natoms += (molecules[imolty].ncount)*(molecules[imolty].nunit);
    }
    if(myid==0)
        printf("Number of atoms from input parameter file:%d\n",sys_pars->natoms);
}

void Input::print_message(int read_result, const char* sec_name)
{
    switch(read_result)
    {
        case 0:
            printf("Finish reading sec %s\n", sec_name);
            break;
        case 1:
            printf("Didnt  find  sec %s\n", sec_name);
            break;
        case 2:
            printf("Missing matching end for %s\n", sec_name);
   }
}
//0: read success 1:didnt find section 2:missing END
int Input::read_section(const char* sec_name)
{
    int find_sec=0, end_read=0, myid = frame->myid;
    char buf[MAXLEN], tmp[MAXLEN];

    if(myid==0)
    {
        rewind(inpara);
        while(fgets(buf,MAXLEN,inpara)!=NULL)
        {
            strip_newline(buf,tmp);
            if(strcmp(tmp,sec_name)==0)
            {
                find_sec = 1;
                break;
            }
        }
    }

    MPI_Bcast(&find_sec,1,MPI_INT,0,MPI_COMM_WORLD);

    if(!find_sec) return 1;

    if(strcmp(sec_name,"FILE")==0)
    {
        file_pars = new File_para;
    }
    else if(strcmp(sec_name,"BUBBLE")==0)
    {
        bubble_pars = new Bubble_para;
    }
    else if(strcmp(sec_name,"SHELL")==0)
    {
        shell_pars = new Shell_para;
    }
    else if(strcmp(sec_name,"FIELD")==0)
    {
        field_pars = new Field_para;
    }
    else if(strcmp(sec_name,"SYSTEM")==0)
    {
        sys_pars = new Sys_para;
    }
    else
    {
        if(myid==0)
             printf("Unknown sec name in read %s!\n",sec_name);
        exit(1);
    }

    memset(buf,0,MAXLEN);

    int m;                                    //number of characters read
    if(myid==0)
    {
        m = 0;
        char end_string[32];
        strcpy(end_string, "END_");
        strcat(end_string, sec_name);
        while(fgets(&buf[m],MAXLEN,inpara)!=NULL)
        {
            if(strstr(&buf[m],"#")!= NULL) continue; //skip comment
            strip_newline(&buf[m],tmp);
            if(strcmp(tmp,end_string)==0)
            {
                end_read = 1;
                break;
            }   
            int len = strlen(&buf[m]);
            m  += len;
        }
        if(end_read)
        {
            buf[m] = '\0';
        }
    }
    
    MPI_Bcast(&end_read,1,MPI_INT,0,MPI_COMM_WORLD); 
    MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD); 
    MPI_Bcast(buf,m,MPI_CHAR,0,MPI_COMM_WORLD);

    if(end_read)
    {
        set_sec_para(sec_name,buf);
    }
    else
    {
        return 2;
    }

    return 0;
}

void Input::set_sec_para(const char* sec_name,  char* lines)
{

    char* current, *next, *pch, *end;
    char substr[256];
    int len = strlen(lines), sublen;

    current = lines;
    end = &lines[len-1];

    //break into sentences
    while(current<=end) 
    {
        next = strchr(current,'\n');
        sublen = strlen(current) - strlen(next);   
        memcpy(substr,current,sublen);
        substr[sublen] = '\0';
        //brean into tokens, max 16 tokens, each of maxsize 32
        char tokens[16][32];
        int ntokens = 0;
        pch = strtok(substr," ");
        while(pch!=NULL)
        {
            strcpy(tokens[ntokens],pch);
            ntokens++;
            pch = strtok(NULL, " ");
        }
        if(ntokens>0) set_line_para(sec_name,tokens, ntokens);
        current = next + 1;
    }
}

void Input::set_line_para(const char* sec_name, char tokens [16][32], int ntokens)
{  
    int myid = frame->myid;
    const char* first = tokens[0];

    if(strcmp(sec_name,"FILE")==0)
    {
        if(strcmp(first,"traj_file")==0)
        {
            assert(ntokens==2);
            strcpy(file_pars->traj_file,tokens[1]);
        }
        else if(strcmp(first,"maxframe")==0)
        {
            assert(ntokens==2);
            file_pars->maxframe = atoi(tokens[1]);
        }
        else
        {
            if(myid==0) 
                printf("Unknown argument %s\n",first);
            exit(1);
        }
    }
    else if(strcmp(sec_name,"BUBBLE")==0)
    {
        if(strcmp(first,"lbubble")==0)
        {
            assert(ntokens==4);
            bubble_pars->lbubble[0] = atoi(tokens[1]); 
            bubble_pars->lbubble[1] = atoi(tokens[2]);
            bubble_pars->lbubble[2] = atoi(tokens[3]);
        }
        else if(strcmp(first,"bead_type")==0)
        {
            for(int i=1;i<ntokens;i++)
            {
                (bubble_pars->bead_type).insert(atoi(tokens[i]));
            }
        }
        else if(strcmp(first,"bubble_mesh")==0)
        {
            assert(ntokens==4);
            bubble_pars->bubble_mesh[0] = atoi(tokens[1]);    
            bubble_pars->bubble_mesh[1] = atoi(tokens[2]);
            bubble_pars->bubble_mesh[2] = atoi(tokens[3]);
        }
        else if(strcmp(first,"ncut")==0)
        {
            assert(ntokens==2);
            bubble_pars->ncut = atoi(tokens[1]);
        }
        else if(strcmp(first,"nliqcut")==0)
        {
            assert(ntokens==2);
            bubble_pars->nliqcut = atoi(tokens[1]);
        }
        else if(strcmp(first,"rcut")==0)
        {
            assert(ntokens==2);
            bubble_pars->rcut = atof(tokens[1]);
        }
        else 
        {
            if(myid==0) 
                printf("Unknown argument %s\n",first);
            exit(1);
        }
    }
    else if(strcmp(sec_name,"SHELL")==0)
    {
        if(strcmp(first,"lshell")==0)
        {
            assert(ntokens==4);
            shell_pars->lshell[0] = atoi(tokens[1]);
            shell_pars->lshell[1] = atoi(tokens[2]);
            shell_pars->lshell[2] = atoi(tokens[3]);
        }
        else if(strcmp(first,"del")==0)
        {   
            assert(ntokens==2);
            shell_pars->del = atoi(tokens[1]);
        }
        else if(strcmp(first,"nshell")==0)
        {
            assert(ntokens==2);
            shell_pars->nshell= atoi(tokens[1]);
        }
        else if(strcmp(first,"delr")==0)
        {
            assert(ntokens==2);
            shell_pars->delr = atof(tokens[1]);
        }
        else if(strcmp(first,"delv")==0)
        {
            assert(ntokens==2);
            shell_pars->delv = atof(tokens[1]);
        }
        else if(strcmp(first,"deln")==0)
        {
            assert(ntokens==2);
            shell_pars->deln = atoi(tokens[1]);
        }
        else 
        {
            if(myid==0) 
                printf("Unknown argument %s\n",first);
            exit(1);
        }
    }
    else if(strcmp(sec_name,"FIELD")==0)
    {
        if(strcmp(first,"lfield")==0)
        {
            assert(ntokens==4);
            field_pars->lfield[0] = atoi(tokens[1]);
            field_pars->lfield[1] = atoi(tokens[2]);
            field_pars->lfield[2] = atoi(tokens[3]);
        }
        else if(strcmp(first,"field_mesh")==0)
        {
            assert(ntokens==4);
            field_pars->field_mesh[0] = atoi(tokens[1]);
            field_pars->field_mesh[1] = atoi(tokens[2]);
            field_pars->field_mesh[2] = atoi(tokens[3]);
        }
        else
        {
            if(myid==0) 
                printf("Unknown argument %s\n",first);
            exit(1);
        }
    }
    else if(strcmp(sec_name,"SYSTEM")==0)
    {
        if((strcmp(first,"boxl")==0))
        {
            assert(ntokens==4);
            sys_pars->boxl[0] = atof(tokens[1]);
            sys_pars->boxl[1] = atof(tokens[2]);
            sys_pars->boxl[2] = atof(tokens[3]);
        }
        else if(strcmp(first,"nmolty")==0)
        {
            assert(ntokens==2);
            sys_pars->nmolty = atoi(tokens[1]);
            sys_pars->init();
            sys_pars->init_flag = 1; 
        }
        else if(strcmp(first,"molty")==0)
        { 
            assert(ntokens>=8);
            assert(sys_pars->init_flag==1);
            Molecule* molecule_array = sys_pars->all_molecule;
            int nread = sys_pars->nmolty_read;
            molecule_array[nread].nunit = atoi(tokens[1]);
            molecule_array[nread].init();
            molecule_array[nread].init_flag = 1;
            molecule_array[nread].ncount = atoi(tokens[2]);
            molecule_array[nread].fix_dof = atoi(tokens[3]);
            int nunit = molecule_array[nread].nunit;
            int start = 4;
            for(int i=0; i<nunit; i++) {molecule_array[nread].beadmass[i]= atof(tokens[start+i]);}   
            start += nunit;             
            for(int i=0; i<nunit; i++) {molecule_array[nread].atomtype[i]= atoi(tokens[start+i]);}
            start += nunit;
            for(int i=0; i<nunit; i++) {molecule_array[nread].chemid[i]= atoi(tokens[start+i]);}
            start += nunit;
            for(int i=0; i<nunit; i++) {strcpy(molecule_array[nread].atomname[i],tokens[start+i]);} 
            sys_pars->nmolty_read += 1;
        }
        else
        {
            if(myid==0) 
                printf("Unknown argument %s\n",first);
            exit(1);
        }
    }
    else
    {
        if(myid==0)
             printf("Unknown sec name in read %s!\n",sec_name);
        exit(1);
    }
}

void Input::strip_newline(const char* oldline, char* newline)
{
    int len = strlen(oldline); 
    strcpy(newline,oldline);
    if( len>0 && newline[len-1]=='\n')
            newline[len-1] = '\0';
}

void Input::print_tokens(char tokens[16][32], int ntokens)
{
    printf("ntokens:%d\n",ntokens);
    for(int i=0;i<ntokens;i++)
    {
        printf("%s \n",tokens[i]);
    }
}
