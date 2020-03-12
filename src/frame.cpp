#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <queue>
#include <vector>
#include <mpi.h>

#include "frame.h"
#include "domain.h"
#include "input.h"
#include "atom.h"
#include "bubble.h"
#include "function.h"

#define MAXLEN  1024         //used in reading small sections
#define MAXLINE 256          //used in reading atom info 
#define CHUNK  1024          //used in reading atom info
#define NWORDS  11           //number of columns in atom

Frame::Frame(const char* iinpara_file)
{
    intraj = NULL;
    input = NULL;
    domain = NULL;
    atom = NULL;
    bubble = NULL;
    shell = NULL;
    strcpy(inpara_file,iinpara_file);
}

Frame::~Frame()
{
    if(input)
    {
        delete input;
        input = NULL;
    }

    if(domain)
    {
         delete domain;
         domain = NULL;
    }

    if(atom) 
    {
        delete atom;
        atom = NULL;
    } 

    if(bubble)
    {
        delete bubble;
        bubble = NULL;
    }

    if(shell)
    {
        delete shell;
        shell = NULL;
    }

    if(intraj)     fclose(intraj);
    if(cart)       MPI_Comm_free(&cart);
}

//set up proc grid
//allocate memory for pointers
void Frame:: init()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    input = new Input(this);     
    input->init();               
    if(myid==0)
    {   
        printf("Finish reading input parameter file\n");
    }

    set_procgrid();   //need boxl in input file to set procgrid

    if(myid==0)
    {
        printf("Topology of MPI tasks:%d %d %d\n", procgrid[0], procgrid[1],procgrid[2]);
    }

    int periods[3], reorder=1;
    periods[0] = periods[1] = periods[2] = 1; 

    MPI_Cart_create(MPI_COMM_WORLD,3,procgrid,periods,reorder,&cart);
    MPI_Cart_get(cart,3,procgrid,periods,mycor);
    MPI_Cart_rank(cart,mycor,&myid_3D);

    if(myid==0)
    {    
        intraj = fopen(input->file_pars->traj_file,"r");
        if(!intraj)
        {   
            printf("Failed to open traj file\n");
            exit(1);
        }
    }

    domain = new Domain(this);
    atom = new Atom(this);
    bubble = new Bubble(this);
    shell = new Shell(this);
    field = new Field(this);
}

void Frame:: set_procgrid()
{
    int** factors;
    int npossible = factor(nproc,NULL);
    allocate_2D<int>(factors,npossible,3,"factors");
    npossible = factor(nproc,factors);
    if(npossible==0)
    {
        if(myid==0)
            printf("Bad choice of nproc!\n");
        exit(1);
    }
    best_factors(npossible,factors,procgrid);
    free_2D<int> (factors);
}

//read head part of one frame
//all information in domain is updated
int Frame:: read_head()
{
    int read_error = 0, m;
    char buf[MAXLEN];
    if(myid==0)
    {
        m = 0;
        int nline = 0;
        while(fgets(&buf[m],MAXLEN,intraj)!=NULL)
        {
            nline++;
            int len = strlen(&buf[m]);
            m += len;
            if(nline==9)
                break;
        }
        if(nline==9)
        {
            buf[m] ='\0';
        }
        else
        {
            read_error = 1;
            printf("Error in reading head! \n");
        }
    }

    MPI_Bcast(&read_error,1,MPI_INT,0,MPI_COMM_WORLD);
    if(read_error)
        return 1;
    MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(buf,m,MPI_CHAR,0,MPI_COMM_WORLD);

    char substr[256];
    char* start=buf, *end;
    int sublen; 
    end = strchr(start,'\n'); 

    //line 2: timestep
    start = end + 1;
    end = strchr(start,'\n');  
    sublen = strlen(start) - strlen(end);
    memcpy(substr,start,sublen);
    substr[sublen] = '\0';
    index = atoi(substr);

    //line 4: natoms
    for(int i=0;i<2;i++)
    {
        start = end + 1;
        end = strchr(start,'\n');
    }
    sublen = strlen(start) - strlen(end);
    memcpy(substr,start,sublen);
    substr[sublen] = '\0';
    int natom_traj = atoi(substr); 

    if(natom_traj!=input->sys_pars->natoms)
    {
        if(myid==0)
        {
            printf("Number of atoms not consistent!\n");
        }
        return 1;
    }

    //line 6,7,8 box boundary 
    start = end + 1;
    end = strchr(start,'\n');
    for(int i=0;i<3;i++)
    {
        start = end + 1;
        end = strchr(start,'\n');
        sublen = strlen(start) - strlen(end);
        memcpy(substr,start,sublen);
        substr[sublen] = '\0';
        char* pch = strtok(substr," ");
        domain->boxlo[i] = atof(pch);
        pch = strtok(NULL, " ");
        domain->boxhi[i] = atof(pch);
        domain->boxlen[i] = domain->boxhi[i] - domain->boxlo[i];
     }

     //set sublo & subhi & sublen
     for(int i=0;i<3;i++)
     {
         domain->sublen[i] = domain->boxlen[i]/procgrid[i];
         domain->sublo[i] = domain->boxlo[i] + mycor[i]*domain->sublen[i];
         domain->subhi[i] = domain->sublo[i] + domain->sublen[i];
     }
     return 0;
}

//read atom information from one frame
int Frame:: read_atom()
{
    char buf[MAXLINE*CHUNK];
    int nchunk, nread, read_error = 0;

    atom->init();           //allocate memory for local atoms

    nread = 0;
    int natoms = atom->natoms;
    while(nread<natoms)
    {    
        nchunk = std::min(natoms-nread,CHUNK);
        int m;
        if(myid==0)
        {
            m = 0;
            for(int i=0;i<nchunk;i++)
            {
                if(!fgets(&buf[m],MAXLINE,intraj))
                {
                    printf("Abnormal end in reading atoms\n");
                    read_error = 1;
                    break;
                }
                m += strlen(&buf[m]);
            }
            if(m)
            {
                buf[m] = '\0';
            }
        }

        MPI_Bcast(&read_error,1,MPI_INT,0,MPI_COMM_WORLD);
        if(read_error)
            return 1;
        MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(buf,m,MPI_CHAR,0,MPI_COMM_WORLD);
        set_atom_info(nchunk,buf);
        nread += nchunk;
    }

    int n = atom->nlocal;
    int sum;
    MPI_Allreduce(&n,&sum,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(sum!=atom->natoms)
    {
        if(myid==0)
        {
            printf("Local number of atoms doenst adds up to correct number! %d != %d\n",sum,atom->natoms);
        }
        return 1;
    }
    return 0;
}

void Frame::set_atom_info(int nline, char* buf)
{
    int atype, aid;
    char* next;
    char** tokens = new char*[NWORDS];

    double recoord[3],coord[3],vel[3],force[3],procgridsize[3];
    int procid[3];

    double* boxlo = domain->boxlo;
    double* boxhi = domain->boxhi;

    double boxlx = domain->boxlen[0];
    double boxly = domain->boxlen[1];
    double boxlz = domain->boxlen[2];

    for(int i=0;i<3;i++)
    {
        procgridsize[i] = 1.0/procgrid[i];
    }

    for(int i=0;i<nline;i++)
    {   
        next = strchr(buf,'\n');
        tokens[0] = strtok(buf," ");
        for(int m=1;m<NWORDS;m++)
        {
            tokens[m] = strtok(NULL," ");  
        }

        aid = atoi(tokens[0]);
        atype = atoi(tokens[1]);

        recoord[0] = atof(tokens[2]);
        recoord[1] = atof(tokens[3]);
        recoord[2] = atof(tokens[4]);

        double lo[3] = {0.,0.,0.};
        double hi[3] = {1.,1.,1.};
        check_bound(recoord,lo,hi);                //map every atom within box boundary

        procid[0] = static_cast<int>(recoord[0]/procgridsize[0]);
        procid[1] = static_cast<int>(recoord[1]/procgridsize[1]);
        procid[2] = static_cast<int>(recoord[2]/procgridsize[2]);

        id_pbc(procid[0],procgrid[0]);
        id_pbc(procid[1],procgrid[1]);
        id_pbc(procid[2],procgrid[2]);
        
        //get my atoms
        if (mycor[0]==procid[0]&&mycor[1]==procid[1]&&mycor[2]==procid[2]) 
        {
            //convert to absolute unit
            coord[0] = boxlx*recoord[0] + boxlo[0];
            coord[1] = boxly*recoord[1] + boxlo[1];
            coord[2] = boxlz*recoord[2] + boxlo[2];
            vel[0] = atof(tokens[5]);
            vel[1] = atof(tokens[6]);
            vel[2] = atof(tokens[7]);
            force[0] = atof(tokens[8]);
            force[1] = atof(tokens[9]);
            force[2] = atof(tokens[10]);
            atom->add_atom_local(aid, atype, coord, vel, force);
        }
        buf = next+1;
    }

    delete [] tokens;
}

int Frame:: factor(int n, int** factors)
{
  int i,j,k,nyz,bubble,bnx,bny,bnz,field,fnx,fny,fnz;
  int* lbubble, *lfield;

  lbubble = input->bubble_pars->lbubble;
  lfield = input->field_pars->lfield;

  bubble = (lbubble[0] || lbubble[1] || lbubble[2]);
  field = (lfield[0] || lfield[1] || lfield[2]);
 
  if(bubble)
  {
      bnx = input->bubble_pars->bubble_mesh[0];
      bny = input->bubble_pars->bubble_mesh[1];
      bnz = input->bubble_pars->bubble_mesh[2];
  }

  if(field)
  {
      fnx = input->field_pars->field_mesh[0];
      fny = input->field_pars->field_mesh[1];
      fnz = input->field_pars->field_mesh[2];
  }

  int m = 0;
  for (i = 1; i <= n; i++) {
    if (n % i || (bubble && bnx%i) || (field && fnx%i)) continue;
    nyz = n/i;
    for (j = 1; j <= nyz; j++) {
      if (nyz % j || (bubble && bny%j) || (field && fny%j)) continue;
      k = nyz/j;
      if((bubble && bnz%k) || (field && fnz%j)) continue;
      if (factors) {
        factors[m][0] = i;
        factors[m][1] = j;
        factors[m][2] = k;
      }
      m++;
    }
  }

  return m;
}

int Frame:: best_factors(int npossible, int **factors, int *best)
{
  int index;
  double surf, boxlx,boxly,boxlz, area[3], bestsurf;

  boxlx = input->sys_pars->boxl[0];
  boxly = input->sys_pars->boxl[1];
  boxlz = input->sys_pars->boxl[2];
  area[0] = boxlx*boxly;
  area[1] = boxlx*boxlz;
  area[2] = boxly*boxlz;

  bestsurf = 2.0 * (area[0]+area[1]+area[2]);

  for (int m = 0; m < npossible; m++) {
    surf = area[0]/factors[m][0]/factors[m][1] +
      area[1]/factors[m][0]/factors[m][2] +
      area[2]/factors[m][1]/factors[m][2];
    if (surf < bestsurf) {
      bestsurf = surf;
      best[0] = factors[m][0];
      best[1] = factors[m][1];
      best[2] = factors[m][2];
      index = m;
    }
  }

  return index;
}
