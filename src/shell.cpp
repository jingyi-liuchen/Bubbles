#include <cmath>
#include <mpi.h>
#include "frame.h"
#include "shell.h"
#include "function.h"
#include "constants.h"
#include "molecule.h"

#define VOLCUT 10      //lower bound of maxvol for setting cms for shell
#define RHOLIQ 1000.0   //bulk liquid water density
#define MAXDEV 0.05    //maximum deviation from deln 

Shell::Shell(class Frame* iframe)
{
    shellr = NULL;
    shellT = NULL;
    shellvr = NULL;
    shellnum = NULL;

    frame = iframe;  
    nshell = frame->input->shell_pars->nshell;
    del = frame->input->shell_pars->del;
}

Shell::~Shell()
{
    if(shellr)   free_1D<double>(shellr);
    if(shellT)   free_1D<double>(shellT);
    if(shellvr)  free_1D<double>(shellvr);
    if(shellnum) free_2D<int>(shellnum);    
}

//driver function
void Shell::cal_shell(FILE* outshell)
{
    int* lshell = frame->input->shell_pars->lshell;
    if(!lshell[0] && !lshell[1] && !lshell[2]) return;

    int myid = frame->myid;
    int maxvol = frame->bubble->maxvol;
    double* maxcms = frame->bubble->maxcms;

    if(maxvol < VOLCUT)
    {
        cms[0] = cms[1] = cms[2] = 0.5;
    }
    else
    {
       cms[0] = maxcms[0];
       cms[1] = maxcms[1];
       cms[2] = maxcms[2];
    }

    //convert to absolute units
    cms[0] *= frame->domain->boxlen[0];
    cms[1] *= frame->domain->boxlen[1];
    cms[2] *= frame->domain->boxlen[2];

    cal_shellr();
    cal_shellprop();
    if(myid==0) output_shell(outshell);
    cleanup_shell();
}

void::Shell::cal_shellprop()
{
    int myid = frame->myid;
    int* lshell = frame->input->shell_pars->lshell;
    int nmolty = frame->input->sys_pars->nmolty;
    int nlocal = frame->atom->nlocal;
    double** x = frame->atom->x;
    double** v = frame->atom->v;
    class Molecule* all_molecule = frame->input->sys_pars->all_molecule;
    int* type = frame->atom->type;

    int** shellnum_local;
    double **shellvcm_local, **shellvcm, *shellmass_local, *shellmass, *shellT_local, *shellvr_local;
    int* shellid;

    allocate_1D<int>(shellid,nlocal,"shellid");
    allocate_2D<int>(shellnum,nshell,nmolty,"shellnum");
    allocate_2D<int>(shellnum_local,nshell,nmolty,"shellnum_local");

    for(int i=0;i<nmolty;i++)
    {
        for(int j=0;j<nshell;j++)
        {
            shellnum_local[i][j] = 0;
        } 
    }

    if(lshell[1])  //shell temperature 
    {
        allocate_1D<double>(shellT,nshell,"shellT");
        allocate_1D<double>(shellmass,nshell,"shellmass");
        allocate_2D<double>(shellvcm,nshell,3,"shellvcm");

        allocate_1D<double>(shellT_local,nshell,"shellT_local");
        allocate_1D<double>(shellmass_local,nshell,"shellmass_local");
        allocate_2D<double>(shellvcm_local,nshell,3,"shellvcm_local");

        for(int i=0;i<nshell;i++)
        {
            shellmass_local[i] = 0.;
            shellT_local[i] = 0.;
            for(int j=0;j<3;j++)
            {
                shellvcm_local[i][j] = 0.;
            }
        }
    }

    if(lshell[2])  //shell radial velocity
    {
        if(!lshell[1]) 
        {
            allocate_1D<double>(shellmass,nshell,"shellmass");
            allocate_1D<double>(shellmass_local,nshell,"shellmass_local");
        }
        allocate_1D<double>(shellvr,nshell,"shellvr");  
        allocate_1D<double>(shellvr_local,nshell,"shellvr_local");

        for(int i=0;i<nshell;i++)
        {
            shellvr_local[i] = 0.;
            shellmass_local[i] = 0.;
        }
    }

    //calculate shell properties
    double* boxlen = frame->domain->boxlen;
    for(int i=0;i<nlocal;i++)
    {
        double distsq = cal_distsq(x[i],cms,boxlen);
        double dist = sqrt(distsq);
        int ishellid  = cal_shellid(dist);   
        shellid[i] = ishellid;
        if(ishellid < nshell)
        {
            int molty = type[i];
            shellnum_local[ishellid][molty-1]++;
            double beadmass = all_molecule[molty-1].beadmass[0];
            if(lshell[1])
            {
                shellvcm_local[ishellid][0] += beadmass * v[i][0];
                shellvcm_local[ishellid][1] += beadmass * v[i][1];
                shellvcm_local[ishellid][2] += beadmass * v[i][2];
                shellmass_local[ishellid] += beadmass;
            }
            if(lshell[2])
            {
                double delx[3], vr;
                delx[0] = x[i][0] - cms[0];
                delx[1] = x[i][1] - cms[1];
                delx[2] = x[i][2] - cms[2];
                mimage(delx,boxlen);
                vr = v[i][0] * delx[0] + v[i][1] * delx[1] + v[i][2] * delx[2];
                vr /= dist;
                shellvr_local[ishellid] += beadmass * vr; 
            }
        }
    }

    MPI_Allreduce(shellnum_local[0],shellnum[0],nmolty*nshell,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    //calculate shell temperature
    if(lshell[1])
    {
         MPI_Allreduce(shellmass_local,shellmass,nshell,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         MPI_Allreduce(shellvcm_local[0],shellvcm[0],nshell*3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
         for(int i=0;i<nshell;i++)
         {
             if(shellmass[i]>0)
             {
                 for(int j=0;j<3;j++)
                 {
                     shellvcm[i][j] /= shellmass[i];
                 }
             }
         }

         for(int i=0;i<nlocal;i++)
         {
             int ishellid = shellid[i];
             if(ishellid < nshell)
             {
                 double vx_sub, vy_sub, vz_sub;
                 int molty = type[i];                
                 double beadmass = all_molecule[molty-1].beadmass[0];
                 vx_sub = v[i][0] - shellvcm[ishellid][0];
                 vy_sub = v[i][1] - shellvcm[ishellid][1];
                 vz_sub = v[i][2] - shellvcm[ishellid][2]; 
                 shellT_local[ishellid] += beadmass * (vx_sub*vx_sub + vy_sub*vy_sub + vz_sub*vz_sub);
             }
         }

         MPI_Allreduce(shellT_local,shellT,nshell,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

         for(int i=0;i<nshell;i++)
         {
             int dof = 0;
             for(int j=0;j<nmolty;j++)
             {
                 dof += shellnum[i][j] * 3;
             }
             if(dof > 0)
             {
                 double tfactor = constants::mvv2e /(dof*constants::boltz);
                 shellT[i] *= tfactor;                
             }
             else
             {
                 shellT[i] = -1.;
             }
         }
    }

    if(lshell[2])
    {
        if(!lshell[1])
        {
            MPI_Allreduce(shellmass_local,shellmass,nshell,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        }   
        MPI_Allreduce(shellvr_local,shellvr,nshell,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

        for(int i=0;i<nshell;i++)
        {
            if(shellmass[i] > 0)
            {
                shellvr[i] = shellvr[i]/shellmass[i]*1E5; //convert to m/s
            }
        }
    }

    free_2D<int>(shellnum_local);
    free_1D<double>(shellT_local);
    free_1D<double>(shellmass_local);
    free_1D<double>(shellmass);
    free_2D<double>(shellvcm);
    free_2D<double>(shellvcm_local);
    free_1D<double>(shellvr_local);
    free_1D<int>(shellid);
}

int Shell:: cal_shellid(double dist)
{
    double distsq = dist * dist;
    double rtov =  4.0/3*M_PI;

    int shellid;
    if(del==0)
    {
        double delr = frame->input->shell_pars->delr;
        shellid = static_cast<int>(dist/delr);
    }
    else if(del==1)
    {
        double delv = frame->input->shell_pars->delv;
        double v = 4.0/3*M_PI * distsq * dist;
        shellid = static_cast<int>(v/delv);
    }
    else //dn = constant
    {
        double roff = shellr[n_inner_shell-1];
        if(dist <= roff)
        {
            if(n_inner_shell==1)
            {
                shellid = 0;
            }
            else
            {
                for(int i=1;i<n_inner_shell;i++)
                {
                    if(dist <= shellr[i] && dist > shellr[i-1])
                    {
                        shellid = i;
                        break;
                    }
                }
            }
        }        
        else
        {
            int deln = frame->input->shell_pars->deln;
            double delv = deln/constants::avogadro*0.018/RHOLIQ*1E30;                     
            double roff3 = roff * roff * roff;
            double r3 = distsq * dist;
            double vdiff = rtov * (r3 -roff3);
            shellid = static_cast<int>(vdiff/delv) + n_inner_shell;
        }
   }

   return shellid;
}

void Shell:: cal_shellr()
{
     allocate_1D<double>(shellr,nshell,"shellr");
     double vtor = 3./4/M_PI;

     if(del==0) //dr = constant
     {
         double delr = frame->input->shell_pars->delr;    
         for(int i=0;i<nshell;i++)
         {
             shellr[i] = (i+1)*delr;
         }
     }
     else if(del==1) // dV = constant
     {
         double delv = frame->input->shell_pars->delv;
         double vsum = delv;
         for(int i=0;i<nshell;i++)
         {
              shellr[i] = pow(vsum*vtor,1.0/3);
              vsum += delv;
         }
     }
     else if(del==2) //dN = constant
     {
         cal_shellr_deln();
     }
     else
     {
         printf("Unknown value for del in shell pars\n");
         exit(1);
     }
}

void Shell::cal_shellr_deln()
{
    int myid = frame->myid;
    int maxvol = frame->bubble->maxvol; 
    int deln = frame->input->shell_pars->deln;
    int nlocal = frame->atom->nlocal;
    double **x = frame->atom->x, vtor = 3./4/M_PI, rtov=1.0/vtor, third=1.0/3;

    double meshvol = 1.;
    for(int i=0;i<3;i++)
    {
        double meshlen = frame->domain->boxlen[i]/frame->input->bubble_pars->bubble_mesh[i];
        meshvol *= meshlen;
    }

    double v0 = maxvol * meshvol;                                   //current bubble volume
    double delv = deln/constants::avogadro*0.018/RHOLIQ*1E30;       //increment of v assuming bulk density 

    if( maxvol < 10)                                         //close to bulk situation
    {   
        double vsum = v0;
        for(int i=0;i<nshell;i++)
        {   
            vsum += delv;
            shellr[i] = pow(vsum*vtor,third);
        }
        n_inner_shell = 1;
        return;
    }

    double r0 = pow(v0*vtor,third);          //current bubble radius
    double rminsq_local = r0*r0;
    double rmax = r0 + 5.0;                  //max radius for high-resolution histogram 
    double rmaxsq = rmax * rmax;
    double rminsq;
    int ninner_local=0,ninner;

    std::vector<double> inner_dist_local;
    for(int i=0;i<nlocal;i++)
    {
        double distsq = cal_distsq(x[i],cms,frame->domain->boxlen);
        if(distsq < rmaxsq)
        {   
            inner_dist_local.push_back(sqrt(distsq));
            ninner_local++;
            if(distsq < rminsq_local)
            {   
                rminsq_local = distsq;
            }
        }
    }

    MPI_Allreduce(&ninner_local,&ninner,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&rminsq_local,&rminsq,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
 
    double rmin = sqrt(rminsq);  
    shellr[0] = rmin - 0.001;                 //make sure no atom exist when r<rmin

    int low_bound = static_cast<int>((1.0 - MAXDEV)*deln);
    int high_bound = static_cast<int>((1.0 + MAXDEV)*deln);

    if(ninner < low_bound)                    //need to increase rmax
    {   
        int delnn = deln  - ninner;
        double delvv = delnn/constants::avogadro*0.018/RHOLIQ*1E30;
        double vmax = rtov*rmaxsq*rmax + delvv;
        rmax = pow(vmax*vtor,third);          //new rmax
        shellr[1] = rmax;
        double vsum = vmax;
        for(int i=2;i<nshell;i++)
        {   
            vsum += delv;
            shellr[i] = pow(vsum*vtor,third);
        }
        n_inner_shell = 2;
    }
    else if(ninner > high_bound)              //need further division 
    {
        //communicate to get all inner_dist
        double* inner_dist;
        int* ninner_local_all, *displs;

        allocate_1D<double>(inner_dist, ninner,"inner_dist");
        allocate_1D<int>(ninner_local_all,frame->nproc,"ninner_local_all");
        allocate_1D<int>(displs,frame->nproc,"displs");

        MPI_Allgather(&ninner_local,1,MPI_INT,ninner_local_all,1,MPI_INT,MPI_COMM_WORLD);

        int psum = 0;
        for(int i=0;i<frame->nproc;i++)
        {
            displs[i] = psum ;
            psum += ninner_local_all[i];
        }

        MPI_Allgatherv(inner_dist_local.data(),ninner_local,MPI_DOUBLE,
                      inner_dist,ninner_local_all,displs,MPI_DOUBLE,MPI_COMM_WORLD);

        free_1D<int>(ninner_local_all);
        free_1D<int>(displs);
 
        int nshell_high = ninner/10;
        //int nshell_high = 10000;
        double delr = (rmax-rmin)/nshell_high;
        int* hist;
        allocate_1D<int>(hist,nshell_high,"shellhist");
        for(int i=0;i<nshell_high;i++)
        {
            hist[i] = 0;
        }

        for(int i=0;i<ninner;i++)
        {
            int shellid = static_cast<int>((inner_dist[i]-rmin)/delr);
            if(shellid>=0 && shellid<=nshell_high-1)  hist[shellid]++;
        }

        int period_sum = 0, ncount=0, end_shell, *bounds;
        allocate_1D<int>(bounds,nshell_high,"bounds");
        for(int i=0;i<nshell;i++)
        {
            if(hist[i] > high_bound) break;
            period_sum += hist[i];
            if (period_sum >= low_bound && period_sum <=high_bound)
            {
                bounds[ncount++] = i;
                period_sum = 0;
            }
        }

        if(ncount==0) //something is wrong, didnt find any new boundary
        {
            printf("Error in dividing inner shells! Likely Resolution not high enough!\n");
            exit(1);
        }

        for(int i=1;i<ncount+1;i++)
        {
            shellr[i] = rmin + delr*(bounds[i-1]+1);
        }

        end_shell = ncount;
        n_inner_shell = ncount + 1;

        double r = shellr[end_shell];
        double vsum =  rtov*r*r*r;
        for(int i=end_shell+1;i<nshell;i++)
        {
            vsum += delv;
            shellr[i] = pow(vsum*vtor,third);
        }

        free_1D<int>(hist);
        free_1D<int>(bounds);
        free_1D<double>(inner_dist);
    }
    else //got lucky 
    {
        shellr[1] = rmax;
        double vsum = rtov*rmaxsq*rmax;
        for(int i=2;i<nshell;i++)
        {
            vsum += delv;
            shellr[i] = pow(vsum*vtor,third);
        }
        n_inner_shell = 2;
    }
    return;
}

void Shell::cleanup_shell()
{
    if(shellr)   free_1D<double>(shellr);
    if(shellT)   free_1D<double>(shellT);
    if(shellvr)  free_1D<double>(shellvr);
    if(shellnum) free_2D<int>(shellnum);
}

void Shell::output_shell(FILE* outshell)
{
    int nmolty = frame->input->sys_pars->nmolty;
    fprintf(outshell, "%15.5e ",shellr[0]/2.0);
    if(shellnum)
    {
        for(int j=0;j<nmolty;j++)
        {
            fprintf(outshell, "%15d ",shellnum[0][j]);
        }
    }
    if(shellT)   fprintf(outshell, "%15.5e ",shellT[0]);
    if(shellvr)  fprintf(outshell, "%15.5e\n",shellvr[0]);
 

    for(int i=1;i<nshell;i++)
    {
        fprintf(outshell, "%15.5e ",(shellr[i]+shellr[i-1])/2.0);
        if(shellnum) 
        {
            for(int j=0;j<nmolty;j++)
            {
                fprintf(outshell, "%15d ",shellnum[i][j]);
            }
        }
        if(shellT)   fprintf(outshell, "%15.5e ",shellT[i]);
        if(shellvr)  fprintf(outshell, "%15.5e\n",shellvr[i]);
    }
}
