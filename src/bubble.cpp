#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <queue>
#include <vector>
#include <mpi.h>
#include <Eigen>

#include "frame.h"
#include "function.h"

#define MAXINCELL 5         //max number of atoms inside a cell

using  namespace Eigen;

Bubble:: Bubble(Frame* iframe)
{
    icell = NULL;
    nicell = NULL;
    iucell = NULL;

    allcellty = NULL;
    kappa = NULL;
    ccms = NULL;
 
    frame = iframe;
}

Bubble::~Bubble()
{
    if(icell)      free_1D<int>(icell);
    if(nicell)     free_1D<int>(nicell);
    if(iucell)     free_2D<int>(iucell);

    if(allcellty)  free_1D<int>(allcellty);
    if(kappa)      free_1D<double>(kappa);
    if(ccms)       free_2D<double>(ccms);
}

void Bubble::cal_bubble(FILE* outbubble)
{
    int* lbubble = frame->input->bubble_pars->lbubble;
    if(!lbubble[0] && !lbubble[1] && !lbubble[2]) return;

    int myid = frame->myid;
    set_linkcell();
    cal_cellty();
    if(myid==0)
    {
        BFS();
        int* sortid = sort_vol();
        if(lbubble[1]) cal_cms(sortid);
        if(lbubble[2]) cal_kappa(sortid);
        output_vol(outbubble);
        free_1D<int>(sortid);
    }

    //required for later shell prop calculation
    MPI_Bcast(&maxvol,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(maxcms,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
    cleanup_bubble();
}

void Bubble::set_linkcell()
{
    int* mesh = frame->input->bubble_pars->bubble_mesh;
    double** x = frame->atom->x;
    int nlocal  = frame->atom->nlocal;
    double rcut = frame->input->bubble_pars->rcut;
    double meshsize[3], sublo[3], subhi[3];

    for(int i=0;i<3;i++)
    {
        meshsize[i] = frame->domain->boxlen[i]/mesh[i];
        nlayer[i] = static_cast<int>(rcut/meshsize[i]) + 1;
        sublo[i] = frame->domain->sublo[i];
        subhi[i] = frame->domain->subhi[i];
    }

    //expand the boundary to include nlayer of surface cells
    ncell = 1;
    for(int i=0;i<3;i++)
    {
        ncellx[i] = mesh[i]/frame->procgrid[i] + nlayer[i]*2;
        ncell *= ncellx[i];
        sublo[i] -= meshsize[i]*nlayer[i];
        subhi[i] += meshsize[i]*nlayer[i];
    }

    //allocate space for remote atom coordinates
    frame->atom->init_remote(ncellx,nlayer,frame->procgrid);
    int nremote_max = frame->atom->nremote_max;

    allocate_1D<int>(nicell,ncell,"nicell");
    for(int i=0; i<ncell; i++) {nicell[i]=0;}
    allocate_2D<int>(iucell,ncell,MAXINCELL,"iucell");
    allocate_1D<int>(icell,nlocal+nremote_max,"icell");

    //calculate local atom information
    for(int i=0;i<nlocal;i++)
    {
        int atype = frame->atom->type[i];
        int id_3D[3];
        if(atype != 1) continue;    //only consider water 
        for(int j=0;j<3;j++)
        {
            id_3D[j] = static_cast<int>((x[i][j]-sublo[j])/meshsize[j]);
            id_pbc(id_3D[j],ncellx[j]);
        }
        int id_1D = mapto1D(id_3D, ncellx);
        icell[i] = id_1D;
        iucell[id_1D][nicell[id_1D]] = i;
        nicell[id_1D]++;
        if(nicell[id_1D]>=MAXINCELL)
        {
            printf("Exceed max number of molecules in a cell!\n");
            exit(1);
        }
    }

    comm_surfz();
    comm_surfy();
    comm_surfx();

}

int Bubble::cal_molty(int imol, int& id_1D)
{
    int id_1D_neigh, nlocal = frame->atom->nlocal;
    int nliqcut = frame->input->bubble_pars->nliqcut;
    int id_3D[3],id_3D_neigh[3];
    double* coor = frame->atom->x[imol],*coor_neigh;
    double** x = frame->atom->x, **xre = frame->atom->xre;
    double rcut = frame->input->bubble_pars->rcut;
    double rcutsq = rcut * rcut;

    id_1D = icell[imol];
    mapto3D(id_1D,id_3D,ncellx);
    int nneigh = 0;
    for(int i=-nlayer[0];i<=nlayer[0];i++)
    {
        id_3D_neigh[0] = id_3D[0] + i;
        id_pbc(id_3D_neigh[0],ncellx[0]);
        for(int j=-nlayer[1];j<=nlayer[1];j++)
        {
            id_3D_neigh[1] = id_3D[1] + j;
            id_pbc(id_3D_neigh[1],ncellx[1]);
            for(int k=-nlayer[2];k<=nlayer[2];k++)
            {
                id_3D_neigh[2] = id_3D[2] + k;
                id_pbc(id_3D_neigh[2],ncellx[2]);
                id_1D_neigh = mapto1D(id_3D_neigh,ncellx);
                int nmol = nicell[id_1D_neigh];
                for(int m=0;m<nmol;m++)
                {
                    int jmol = iucell[id_1D_neigh][m];
                    if(imol==jmol) continue;
                    if(jmol<nlocal)
                    {
                        coor_neigh = x[jmol];
                    }
                    else
                    {
                        coor_neigh = xre[jmol-nlocal];
                    }
                    double distsq = cal_distsq(coor, coor_neigh, frame->domain->boxlen);
                    if(distsq < rcutsq)
                    {
                        nneigh++;
                        if(nneigh>=nliqcut)
                        {
                            return 1;
                        }
                    }
                }
            }
        }
    }

    return 0;
}

//cellty 0:vapor 1:liquid
void Bubble::cal_cellty()
{  
    int nliqmol=0, nliqcell=0, nlocal=frame->atom->nlocal, nremote=frame->atom->nremote;
    int nliqcut = frame->input->bubble_pars->nliqcut, ncell = ncellx[0]*ncellx[1]*ncellx[2];
    int* molty,*remolty,*cellty,*cellnnliq;
    int* mesh = frame->input->bubble_pars->bubble_mesh;
    double** x = frame->atom->x, **xre = frame->atom->xre;
    double sublo[3];
    double rcut = frame->input->bubble_pars->rcut;
    double rcutsq = rcut * rcut;
    double meshsize[3];

    for(int i=0;i<3;i++)
    {
        meshsize[i] = frame->domain->boxlen[i]/mesh[i];
        sublo[i] = frame->domain->sublo[i] - meshsize[i]*nlayer[i];
    }

    allocate_1D<int>(molty,nlocal,"molty");
    allocate_1D<int>(remolty,nremote,"remolty");
    allocate_1D<int>(cellty,ncell,"cellty");
    allocate_1D<int>(cellnnliq,ncell,"cellnnliq");

    for(int i=0;i<ncell;i++)
    {
        cellty[i] = 0;
        cellnnliq[i] = 0;
    }

    //mark cells that contain liquid molecule as liquid cells 
    for(int imol=0;imol<nlocal;imol++)
    {
       int atype = frame->atom->type[imol];
       if(atype!=1) continue;
       int id_1D,imolty;
       imolty = cal_molty(imol,id_1D);
       molty[imol] = imolty;
       if(imolty==1)
       {
           nliqmol++;
           if(cellty[id_1D]==0)
           {
               cellty[id_1D] = 1;
               nliqcell++;
           }
       }
    }

    //local liquid molecules
    //mark cells that dont contain liquid molecule but have at least nliqcut
    //nearest neighbor of liquid molecules as liquid cell
    for(int imol=0;imol<nlocal;imol++)
    {
        int atype = frame->atom->type[imol];
        if(atype!=1) continue;        
        if(molty[imol]!=1) continue;
        int id_1D = icell[imol];
        int id_3D[3], id_3D_neigh[3];
        mapto3D(id_1D,id_3D,ncellx);
        for(int i=-nlayer[0];i<=nlayer[0];i++)
        { 
            id_3D_neigh[0] = id_3D[0] + i;
            if(id_3D_neigh[0]<nlayer[0] || id_3D_neigh[0]>ncellx[0]-nlayer[0]-1) 
                continue;
            for(int j=-nlayer[1];j<=nlayer[1];j++)
            {
                id_3D_neigh[1] = id_3D[1] + j;
                if(id_3D_neigh[1]<nlayer[1] || id_3D_neigh[1]>ncellx[1]-nlayer[1]-1) 
                    continue;
                for(int k=-nlayer[2];k<=nlayer[2];k++)
                {
                    id_3D_neigh[2] = id_3D[2] + k;
                    if(id_3D_neigh[2]<nlayer[2] || id_3D_neigh[2]>ncellx[2]-nlayer[2]-1 || (i==0&&j==0&&k==0)) 
                        continue;
                    int id_1D_neigh = mapto1D(id_3D_neigh,ncellx);
                    if(cellty[id_1D_neigh]==0)
                    {
                        double gridcms[3];
                        gridtocms(id_3D_neigh,meshsize,sublo,gridcms);
                        double distsq = cal_distsq(x[imol],gridcms,frame->domain->boxlen);
                        if(distsq<rcutsq)
                        {
                            cellnnliq[id_1D_neigh]++;
                            if(cellnnliq[id_1D_neigh]>=nliqcut)
                            {
                                cellty[id_1D_neigh] = 1;
                                nliqcell++;
                            }
                        }
                    }
                }  
            }
        }
    }

    int n1 = comm_surfz_molty(molty,remolty,0);
    int n2 = comm_surfy_molty(molty,remolty,n1);
    int n3 = comm_surfx_molty(molty,remolty,n2);

    //remote liquid molecules
    for(int imol=0;imol<nremote;imol++)
    {
        if(remolty[imol]!=1) continue;
        int id_1D = icell[imol+nlocal];
        int id_3D[3], id_3D_neigh[3];
        mapto3D(id_1D,id_3D,ncellx);
        for(int i=-nlayer[0];i<=nlayer[0];i++)
        {
            id_3D_neigh[0] = id_3D[0] + i;
            if(id_3D_neigh[0]<nlayer[0] || id_3D_neigh[0]>ncellx[0]-nlayer[0]-1)
                continue;
            for(int j=-nlayer[1];j<=nlayer[1];j++)
            {
                id_3D_neigh[1] = id_3D[1] + j;
                if(id_3D_neigh[1]<nlayer[1] || id_3D_neigh[1]>ncellx[1]-nlayer[1]-1)
                    continue;
                for(int k=-nlayer[2];k<=nlayer[2];k++)
                {
                    id_3D_neigh[2] = id_3D[2] + k;
                    if(id_3D_neigh[2]<nlayer[2] || id_3D_neigh[2]>ncellx[2]-nlayer[2]-1 || (i==0&&j==0&&k==0))
                        continue;
                    int id_1D_neigh = mapto1D(id_3D_neigh,ncellx);
                    if(cellty[id_1D_neigh]==0)
                    {
                        double gridcms[3];
                        gridtocms(id_3D_neigh,meshsize,sublo,gridcms);
                        double distsq = cal_distsq(xre[imol],gridcms,frame->domain->boxlen);
                        if(distsq<rcutsq)
                        {
                            cellnnliq[id_1D_neigh]++;
                            if(cellnnliq[id_1D_neigh]>=nliqcut)
                            {
                                cellty[id_1D_neigh] = 1;
                                nliqcell++;
                            }
                        }
                    }
                }
            }
        }
    }

    free_1D<int>(molty);
    free_1D<int>(remolty);
    free_1D<int>(cellnnliq);

    free_1D<int>(nicell);
    free_2D<int>(iucell);
    free_1D<int>(icell);

    comm_cellty(cellty);

    free_1D<int>(cellty);
}


void Bubble::comm_cellty(const int* cellty)
{
    MPI_Datatype sendarray,recvarray,resized_recvarray;
    int sendstarts[3],sendsubsizes[3],sendsizes[3];
    int recvstarts[3],recvsubsizes[3],recvsizes[3];
    int myid = frame->myid;
    int* mycor = frame->mycor;
    int nproc = frame->nproc;

    if(myid!=0)
    {
        allocate_1D<int>(allcellty,1, "allcellty");
    }
    else
    {
        int ncellall = 1;
        for(int i=0;i<3;i++)
        {   
            ncellall *= frame->input->bubble_pars->bubble_mesh[i];
        }
        allocate_1D<int>(allcellty,ncellall, "allcellty");
    }

    for(int i=0;i<3;i++)
    {
        sendstarts[i] = nlayer[i];
        sendsubsizes[i] = ncellx[i]-2*nlayer[i];
        sendsizes[i] = ncellx[i];
        recvstarts[i] = 0;
        recvsubsizes[i] = sendsubsizes[i];
        recvsizes[i] = frame->input->bubble_pars->bubble_mesh[i];
    }

    MPI_Type_create_subarray(3,sendsizes,sendsubsizes,sendstarts,MPI_ORDER_C,MPI_INT,&sendarray);
    MPI_Type_commit(&sendarray);
    MPI_Type_create_subarray(3,recvsizes,recvsubsizes,recvstarts,MPI_ORDER_C,MPI_INT,&recvarray);
    MPI_Type_commit(&recvarray);
    MPI_Type_create_resized(recvarray,0,sizeof(int),&resized_recvarray);
    MPI_Type_commit(&resized_recvarray);

    int* displs, *recvcounts;
    allocate_1D<int>(displs,nproc,"displs");
    allocate_1D<int>(recvcounts,nproc,"recvcounts");
    for(int i=0;i<nproc;i++)
    {
        recvcounts[i] = 1;
    } 
 
    int mydisp_3D[3],mydisp_1D;
    for(int i=0;i<3;i++)
    {
        int ndisp = frame->input->bubble_pars->bubble_mesh[i]/frame->procgrid[i];
        mydisp_3D[i] = mycor[i]*ndisp;
    } 
    mydisp_1D = mapto1D(mydisp_3D,frame->input->bubble_pars->bubble_mesh);
    MPI_Allgather(&mydisp_1D,1,MPI_INT,displs,1,MPI_INT,MPI_COMM_WORLD);     

    MPI_Gatherv(cellty,1,sendarray,allcellty,recvcounts,displs,resized_recvarray,0,MPI_COMM_WORLD);
    MPI_Type_free(&sendarray); 
    MPI_Type_free(&recvarray);    
    MPI_Type_free(&resized_recvarray);
    free_1D<int>(displs);
    free_1D<int>(recvcounts);
}

void Bubble::BFS()
{
    int* bubble_mesh = frame->input->bubble_pars->bubble_mesh;
    int nmesh = bubble_mesh[0]*bubble_mesh[1]*bubble_mesh[2];
    int* visited;

    vcluster.clear();
    clustermem.clear();
    allocate_1D<int>(visited,nmesh,"visited");
    for(int i=0;i<nmesh;i++) 
    {
        visited[i] = 0;
    }

    ncluster = 0;
    for(int i=0;i<nmesh;i++) 
    {
        if(allcellty[i]==0 && !visited[i])
        {
            int nmem = BFS_visit(i,bubble_mesh,visited);
            vcluster.push_back(nmem);
            ncluster++;
        }
    } 
    free_1D<int>(visited);
}

int Bubble::BFS_visit(int id_1D, const int* mesh, int* visited)
{
    int id_1D_neigh, id_3D[3], id_3D_neigh[3], nmem=0;
    std::queue<int> myq;
    std::vector<int> iclustermem;

    iclustermem.push_back(id_1D);
    myq.push(id_1D);
    visited[id_1D] = 1;
    nmem++;

    while(!myq.empty())
    {
        id_1D = myq.front(); 
        myq.pop();
        mapto3D(id_1D,id_3D,mesh);
        for(int i=-1;i<=1;i++)
        {
            id_3D_neigh[0] = id_3D[0] + i;
            id_pbc(id_3D_neigh[0],mesh[0]);
            for(int j=-1;j<=1;j++)
            { 
                id_3D_neigh[1] = id_3D[1] + j;
                id_pbc(id_3D_neigh[1],mesh[1]);
                for(int k=-1;k<=1;k++)
                {
                    if(i==0&&j==0&&k==0) continue;
                    id_3D_neigh[2] = id_3D[2] + k;
                    id_pbc(id_3D_neigh[2],mesh[2]);

                    id_1D_neigh = mapto1D(id_3D_neigh,mesh);
                    if(!visited[id_1D_neigh] && allcellty[id_1D_neigh]==0)
                    {
                         iclustermem.push_back(id_1D_neigh);
                         myq.push(id_1D_neigh);
                         visited[id_1D_neigh] = 1;
                         nmem++;
                    }
                }
            }
       }
    }
    clustermem.push_back(iclustermem);
    return nmem;
}

int* Bubble::sort_vol()
{
    int ncut = frame->input->bubble_pars->ncut; 
    int vtmp, idtmp;
    int* sortid;
    allocate_1D<int>(sortid,ncluster,"sortid");

    for(int i=0;i<ncluster;i++)
    {
        sortid[i] = i;
    }

    realncut = ncut;
    if(ncluster<ncut)
    {
       realncut = ncluster;
    }

    //sort the first realncut elements
    for(int i=0;i<realncut;i++)
    {
       for(int j=ncluster-1;j>i;j--)
       {
           if (vcluster[j] > vcluster[j-1])
           {
               vtmp = vcluster[j];
               vcluster[j] = vcluster[j-1];
               vcluster[j-1] = vtmp;
               idtmp = sortid[j];
               sortid[j] = sortid[j-1];
               sortid[j-1] = idtmp;
           }
       }
    }

    maxvol = 0;
    if(realncut) maxvol = vcluster[0];
    return sortid;
}

void Bubble::cal_cms(const int* sortid)
{
    if(realncut==0)
    {
        maxcms[0] = maxcms[1] = maxcms[2] = 0.5;
        return;
    }

    double thetax,thetay,thetaz,cosavex,sinavex,cosavey,sinavey,cosavez,sinavez;
    double meshsize[3];
    int* mesh = frame->input->bubble_pars->bubble_mesh;
    double reboxlo[3] = {0.,0.,0.};

    allocate_2D<double>(ccms,realncut,3,"ccms");
    for(int i=0;i<3;i++)
    {
        meshsize[i] = 1.0/mesh[i];
    }

    for(int i=0;i<realncut;i++)
    {   
        int id = sortid[i];
        int nmem = clustermem[id].size();
        cosavex = 0.0;
        sinavex = 0.0;
        cosavey = 0.0;
        sinavey = 0.0;
        cosavez = 0.0;
        sinavez = 0.0;
        for(int j=0;j<nmem;j++)
        {   
            int id_1D = clustermem[id][j];
            int id_3D[3];
            double cms[3];
            mapto3D(id_1D,id_3D,mesh);
            gridtocms(id_3D,meshsize,reboxlo,cms);    //cms in reduced unit
            thetax = cms[0]*2.0*M_PI - M_PI; //circular mean
            thetay = cms[1]*2.0*M_PI - M_PI;
            thetaz = cms[2]*2.0*M_PI - M_PI;
            cosavex += cos(thetax);
            sinavex += sin(thetax);
            cosavey += cos(thetay);
            sinavey += sin(thetay);
            cosavez += cos(thetaz);
            sinavez += sin(thetaz);
        }
        
        cosavex /= nmem;
        sinavex /= nmem;
        cosavey /= nmem;
        sinavey /= nmem;
        cosavez /= nmem;
        sinavez /= nmem;

        ccms[i][0] = (atan2(sinavex,cosavex) + M_PI)/2.0/M_PI;
        ccms[i][1] = (atan2(sinavey,cosavey) + M_PI)/2.0/M_PI;
        ccms[i][2] = (atan2(sinavez,cosavez) + M_PI)/2.0/M_PI;
    }

    maxcms[0] = ccms[0][0];
    maxcms[1] = ccms[0][1];
    maxcms[2] = ccms[0][2];
}

//kappa is calculated after cms
void Bubble:: cal_kappa(const int* sortid)
{
    if(realncut==0)  return;
    
    int* mesh = frame->input->bubble_pars->bubble_mesh;
    double  xzratio = static_cast<double>(mesh[0])/mesh[2]; //correct for using reduced coordinates
    double  yzratio = static_cast<double>(mesh[1])/mesh[2];
    double  meshsize[3], reboxlo[3] = {0.,0.,0.},reboxlen[3] = {1.,1.,1.};

    allocate_1D<double>(kappa,realncut,"kappa");
    for(int i=0;i<3;i++)
    {
        meshsize[i] = 1.0/mesh[i];
    }
    
    for(int i=0;i<realncut;i++)
    {   
        int id = sortid[i];
        int nmem = clustermem[id].size();
        Matrix3d gyraT;
        double lambda1,lambda2,lambda3,lambda_sum;

        for(int l=0;l<3;l++)
        {
            for(int m=0;m<3;m++)
            {
                gyraT(l,m) = 0.0;
            }
        }

        for(int j=0;j<nmem;j++)
        {  
            int id_1D = clustermem[id][j];
            int id_3D[3];
            double cms[3],delx[3];
            mapto3D(id_1D,id_3D,mesh);
            gridtocms(id_3D,meshsize,reboxlo,cms); 
            for(int k=0;k<3;k++)
            {
                delx[k] = cms[k] - ccms[i][k];
            }
            mimage(delx,reboxlen);
            delx[0] *= xzratio;
            delx[1] *= yzratio;
            gyraT(0,0) += delx[0]*delx[0];
            gyraT(1,1) += delx[1]*delx[1];
            gyraT(2,2) += delx[2]*delx[2];
            gyraT(0,1) += delx[0]*delx[1];
            gyraT(0,2) += delx[0]*delx[2];
            gyraT(1,2) += delx[1]*delx[2];
        }
        gyraT(0,0) /= nmem;
        gyraT(1,1) /= nmem;
        gyraT(2,2) /= nmem;
        gyraT(0,1) /= nmem;
        gyraT(0,2) /= nmem;
        gyraT(1,2) /= nmem;

        gyraT(1,0) = gyraT(0,1);
        gyraT(2,0) = gyraT(0,2);
        gyraT(2,1) = gyraT(1,2);

        EigenSolver<MatrixXd> es(gyraT,false);
        lambda1 = (es.eigenvalues())[0].real();
        lambda2 = (es.eigenvalues())[1].real();
        lambda3 = (es.eigenvalues())[2].real();
        lambda_sum = lambda1 + lambda2 + lambda3;
        kappa[i] = 1-3*(lambda1*lambda2+lambda2*lambda3+lambda3*lambda1)/lambda_sum/lambda_sum;
    }
}

void Bubble::output_vol(FILE* outbubble)
{
    fprintf(outbubble,"%d %d\n",frame->index, realncut);
    for(int i=0;i<realncut;i++)
    {
       fprintf(outbubble,"%d\n",vcluster[i]);
       if(ccms)
           fprintf(outbubble,"%10.5e %10.5e %10.5e\n",ccms[i][0],ccms[i][1],ccms[i][2]); 
       if(kappa)
           fprintf(outbubble,"%10.5e\n",kappa[i]);    
    }
}

void Bubble::cleanup_bubble()
{
    if(allcellty)  free_1D<int>(allcellty);
    if(kappa)      free_1D<double>(kappa);
    if(ccms)       free_2D<double>(ccms);
}

int Bubble::comm_surfz_molty(const int* molty, int* remolty, int offset)
{
    int disp = -1, source,dest,nrecv_tot=offset;
    MPI_Status status;
    MPI_Comm cart = frame->cart;
    MPI_Cart_shift(frame->cart,2,disp,&source,&dest);
          
    int nsurfcell = ncellx[0]*ncellx[1]*nlayer[2];
    int nsend_mol=0, nrecv_mol=0;
    int *sendmolty;                          
    int nmaxsend = nsurfcell*MAXINCELL; 

    allocate_1D<int>(sendmolty,nmaxsend, "sendmolty");             
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {
            for(int k = nlayer[2];k<2*nlayer[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    sendmolty[nsend_mol++] = molty[mol_id];
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    offset += nrecv_mol;
    nrecv_tot += nrecv_mol;

    disp = 1;
    MPI_Cart_shift(cart,2,disp,&source,&dest);

    nsend_mol=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {
            for(int k = ncellx[2]-2*nlayer[2];k<ncellx[2]-nlayer[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    sendmolty[nsend_mol++] = molty[mol_id];
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    nrecv_tot += nrecv_mol;
    free_1D<int>(sendmolty);
    return nrecv_tot;
}

int Bubble::comm_surfy_molty(const int* molty, int* remolty, int offset)
{
    int disp = -1, source,dest,nrecv_tot=offset;
    int nlocal = frame->atom->nlocal;
    MPI_Status status;
    MPI_Comm cart = frame->cart;
    MPI_Cart_shift(cart,1,disp,&source,&dest);

    int nsurfcell = ncellx[0]*ncellx[2]*nlayer[1];
    int nsend_mol=0, nrecv_mol=0;
    int *sendmolty;
    int nmaxsend = nsurfcell*MAXINCELL;

    allocate_1D<int>(sendmolty,nmaxsend,"sendmolty");
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {   
        for(int j = nlayer[1];j<2*nlayer[1];j++)
        {   
            for(int k =0;k<ncellx[2];k++)
            {   
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {   
                    int mol_id = iucell[id_1D][m];
                    if(mol_id<nlocal)
                    {
                        sendmolty[nsend_mol++] = molty[mol_id];
                    }
                    else
                    {
                        sendmolty[nsend_mol++] = remolty[mol_id-nlocal];
                    }
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    offset += nrecv_mol;
    nrecv_tot += nrecv_mol;

    disp = 1;
    MPI_Cart_shift(cart,1,disp,&source,&dest);

    nsend_mol=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = ncellx[1]-2*nlayer[1];j<ncellx[1]-nlayer[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    if(mol_id<nlocal)
                    {
                        sendmolty[nsend_mol++] = molty[mol_id];
                    }
                    else
                    {
                        sendmolty[nsend_mol++] = remolty[mol_id-nlocal];
                    }
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    nrecv_tot += nrecv_mol;
    free_1D<int>(sendmolty);
    return nrecv_tot;
}

int Bubble::comm_surfx_molty(const int* molty, int* remolty, int offset)
{
    int disp = -1, source,dest,nrecv_tot=offset;
    int nlocal = frame-> atom->nlocal;
    MPI_Status status;
    MPI_Comm cart = frame->cart;
    MPI_Cart_shift(cart,0,disp,&source,&dest);

    int nsurfcell = ncellx[1]*ncellx[2]*nlayer[0];
    int nsend_mol=0, nrecv_mol=0;
    int *sendmolty;
    int nmaxsend = nsurfcell*MAXINCELL;

    allocate_1D<int>(sendmolty,nmaxsend,"sendmolty");
    for(int i = nlayer[0];i<2*nlayer[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {   
                    int mol_id = iucell[id_1D][m];
                    if(mol_id<nlocal)
                    {   
                        sendmolty[nsend_mol++] = molty[mol_id];
                    }
                    else
                    {   
                        sendmolty[nsend_mol++] = remolty[mol_id-nlocal];
                    }
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    offset += nrecv_mol;
    nrecv_tot += nrecv_mol;

    disp = 1;
    MPI_Cart_shift(cart,0,disp,&source,&dest);

    nsend_mol=0;
    for(int i = ncellx[0]-2*nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    if(mol_id<nlocal)
                    {
                        sendmolty[nsend_mol++] = molty[mol_id];
                    }
                    else
                    {
                        sendmolty[nsend_mol++] = remolty[mol_id-nlocal];
                    }
                }
           }
        }
    }

    MPI_Sendrecv(&nsend_mol,1,MPI_INT,dest,1,&nrecv_mol,1,MPI_INT,source,1,cart,&status);
    MPI_Sendrecv(sendmolty,nsend_mol,MPI_INT,dest,1,&remolty[offset],nrecv_mol,MPI_INT,source,1,cart,&status);

    nrecv_tot += nrecv_mol;
    free_1D<int>(sendmolty);
    return nrecv_tot;
}

void Bubble::comm_surfz()
{
   //first send to -z direction, receive from +z direction    
    int disp = -1, source,dest;
    MPI_Status status;
    MPI_Comm cart = frame->cart; 
    MPI_Cart_shift(cart,2,disp,&source,&dest);
   
    int nsurfcell = ncellx[0]*ncellx[1]*nlayer[2];
    int ncount_cell=0, ncount_mol=0;
    int* sendbuf;                          
    double* sendcoor, *recvcoor;
    double** x = frame->atom->x;
    int nmaxsend = nsurfcell*MAXINCELL; 

    allocate_1D<int>(sendbuf,nsurfcell, "sendbuf");                
    allocate_1D<double>(sendcoor,nmaxsend*3, "sendcoor");            
    allocate_1D<double>(recvcoor,nmaxsend*3, "recvcoor");

    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {
            for(int k = nlayer[2];k<2*nlayer[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);                            
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;   
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    sendcoor[ncount_mol*3] = x[mol_id][0];
                    sendcoor[ncount_mol*3+1] = x[mol_id][1];
                    sendcoor[ncount_mol*3+2] = x[mol_id][2];
                    ncount_mol++;
                }
           }
        }
    }

    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    int nrecv = 0;
    for(int i=0;i<ncount_cell;i++)
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    //process received coordinates: update nicell, icell, iucell, and xre
    int ncount_cell2 = 0, ncount_mol2=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {   
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {   
            for(int k = ncellx[2]-nlayer[2];k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;    
                for(int m=0;m<nmolecell;m++)
                { 
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3];
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2]; 
                    frame->atom->add_atom_remote(coor); 
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

    //send to +z direction, receive from -z direction    
    disp = 1;
    MPI_Cart_shift(cart,2,disp,&source,&dest);
    ncount_cell=0;
    ncount_mol=0;

    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {   
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {   
            for(int k = ncellx[2]-2*nlayer[2];k<ncellx[2]-nlayer[2];k++)
            {   
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;
                for(int m=0;m<nmolcell;m++)
                {   
                    int mol_id = iucell[id_1D][m];
                    sendcoor[ncount_mol*3] = x[mol_id][0];
                    sendcoor[ncount_mol*3+1] = x[mol_id][1];
                    sendcoor[ncount_mol*3+2] = x[mol_id][2];
                    ncount_mol++;
                }
           }
        }
    }

    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    nrecv = 0;
    for(int i=0;i<ncount_cell;i++)
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    ncount_cell2 = 0, ncount_mol2=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = nlayer[1];j<ncellx[1] - nlayer[1];j++)
        {
            for(int k = 0;k < nlayer[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;
                for(int m=0;m<nmolecell;m++)
                {
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3];
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2];
                    frame->atom->add_atom_remote(coor);
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

    free_1D<int>(sendbuf);
    free_1D<double>(sendcoor);
    free_1D<double>(recvcoor);
}

void Bubble::comm_surfy()
{
    //first send to -Y, receive from + Y
    int disp = -1, source,dest;
    MPI_Status status;
    MPI_Comm cart = frame->cart;
    MPI_Cart_shift(cart,1,disp,&source,&dest);

    int nlocal = frame->atom->nlocal;
    int nsurfcell = ncellx[0]*ncellx[2]*nlayer[1];
    int ncount_cell=0, ncount_mol=0;
    int* sendbuf;
    double* sendcoor, *recvcoor;
    double** x = frame->atom->x;
    double** xre = frame->atom->xre;
    int nmaxsend = nsurfcell*MAXINCELL;

    allocate_1D<int>(sendbuf,nsurfcell, "sendbuff");
    allocate_1D<double>(sendcoor,nmaxsend*3, "sendcoor");
    allocate_1D<double>(recvcoor,nmaxsend*3, "recvcoor");

    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = nlayer[1];j<2*nlayer[1];j++)
        {
            for(int k=0;k<ncellx[2];k++) 
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    double* coor;
                    if(mol_id < nlocal)
                    {
                        coor = x[mol_id];
                    }
                    else
                    {
                        coor = xre[mol_id-nlocal];
                    }
                    sendcoor[ncount_mol*3] = coor[0];
                    sendcoor[ncount_mol*3+1] = coor[1];
                    sendcoor[ncount_mol*3+2] = coor[2];
                    ncount_mol++;
                }
           }
        }
    }   
    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    int nrecv = 0;
    for(int i=0;i<ncount_cell;i++) 
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    int ncount_cell2 = 0, ncount_mol2=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {   
        for(int j = ncellx[1]-nlayer[1];j<ncellx[1];j++)
        {
            for(int k = 0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;
                for(int m=0;m<nmolecell;m++)
                {
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3]; 
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2];
                    frame->atom->add_atom_remote(coor);
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

   //second send to +Y, receive from -Y
    disp = 1;
    MPI_Cart_shift(cart,1,disp,&source,&dest);
    ncount_cell=0;
    ncount_mol=0;

    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = ncellx[1]-2*nlayer[1];j<ncellx[1]-nlayer[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    double* coor;
                    if(mol_id < nlocal)
                    {   
                        coor = x[mol_id];
                    }
                    else
                    {   
                        coor = xre[mol_id-nlocal];
                    }
                    sendcoor[ncount_mol*3] = coor[0];
                    sendcoor[ncount_mol*3+1] = coor[1];
                    sendcoor[ncount_mol*3+2] = coor[2];
                    ncount_mol++;
                }
           }
        }
    }
    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    nrecv = 0;
    for(int i=0;i<ncount_cell;i++)
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    ncount_cell2 = 0;
    ncount_mol2=0;
    for(int i = nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = 0;j<nlayer[1];j++)
        {
            for(int k = 0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;
                for(int m=0;m<nmolecell;m++)
                {
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3];
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2];
                    frame->atom->add_atom_remote(coor);
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

    free_1D<int>(sendbuf);
    free_1D<double>(sendcoor);
    free_1D<double>(recvcoor);
}

void Bubble::comm_surfx()
{
    //first send to -X, receive from +X
    int disp = -1, source,dest;
    MPI_Status status;
    MPI_Comm cart = frame->cart;
    MPI_Cart_shift(cart,0,disp,&source,&dest);

    int nlocal = frame->atom->nlocal;
    int nsurfcell = ncellx[1]*ncellx[2]*nlayer[0];
    int ncount_cell=0, ncount_mol=0;
    int* sendbuf;
    double* sendcoor, *recvcoor;
    double** x = frame->atom->x;
    double** xre = frame->atom->xre;
    int nmaxsend = nsurfcell*MAXINCELL;

    allocate_1D<int>(sendbuf,nsurfcell, "sendbuf");
    allocate_1D<double>(sendcoor,nmaxsend*3, "sendcoor");
    allocate_1D<double>(recvcoor,nmaxsend*3, "recvcoor");

    for(int i = nlayer[0];i<2*nlayer[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    double* coor;
                    if(mol_id < nlocal)
                    {
                        coor = x[mol_id];
                    }
                    else
                    {
                        coor = xre[mol_id-nlocal];
                    }
                    sendcoor[ncount_mol*3] = coor[0];
                    sendcoor[ncount_mol*3+1] = coor[1];
                    sendcoor[ncount_mol*3+2] = coor[2];
                    ncount_mol++;
                }
           }
        }
    }
    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    int nrecv = 0;
    for(int i=0;i<ncount_cell;i++)
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    int ncount_cell2 = 0, ncount_mol2=0;
    for(int i = ncellx[0]-nlayer[0];i<ncellx[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k = 0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;
                for(int m=0;m<nmolecell;m++)
                {
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3];
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2];
                    frame->atom->add_atom_remote(coor);
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

    //second send to +X, receive from -X    
    disp = 1;
    MPI_Cart_shift(cart,0,disp,&source,&dest);
    ncount_cell=0;
    ncount_mol=0;

    for(int i = ncellx[0]-2*nlayer[0];i<ncellx[0]-nlayer[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k=0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolcell = nicell[id_1D];
                sendbuf[ncount_cell++] = nmolcell;
                for(int m=0;m<nmolcell;m++)
                {
                    int mol_id = iucell[id_1D][m];
                    double* coor;
                    if(mol_id < nlocal)
                    {
                        coor = x[mol_id];
                    }
                    else
                    {
                        coor = xre[mol_id-nlocal];
                    }
                    sendcoor[ncount_mol*3] = coor[0];
                    sendcoor[ncount_mol*3+1] = coor[1];
                    sendcoor[ncount_mol*3+2] = coor[2];
                    ncount_mol++;
                }
           }
        }
    }
    MPI_Sendrecv_replace(sendbuf,ncount_cell,MPI_INT,dest,1,source,1,cart,&status);
    nrecv = 0;
    for(int i=0;i<ncount_cell;i++)
    {
        nrecv += sendbuf[i];
    }
    MPI_Sendrecv(sendcoor,ncount_mol*3,MPI_DOUBLE,dest,1,recvcoor,nrecv*3,MPI_DOUBLE,source,1,cart,&status);

    ncount_cell2 = 0;
    ncount_mol2=0;
    for(int i =0 ;i<nlayer[0];i++)
    {
        for(int j = 0;j<ncellx[1];j++)
        {
            for(int k = 0;k<ncellx[2];k++)
            {
                int id_3D[3] = {i,j,k};
                int id_1D = mapto1D(id_3D, ncellx);
                int nmolecell = sendbuf[ncount_cell2];
                nicell[id_1D] += nmolecell;
                for(int m=0;m<nmolecell;m++)
                {
                    int atom_index = frame->atom->nlocal + frame->atom->nremote;
                    icell[atom_index] = id_1D;
                    iucell[id_1D][m] = atom_index;
                    double coor[3];
                    coor[0] = recvcoor[ncount_mol2*3];
                    coor[1] = recvcoor[ncount_mol2*3 + 1];
                    coor[2] = recvcoor[ncount_mol2*3 + 2];
                    frame->atom->add_atom_remote(coor);
                    ncount_mol2++;
                }
                ncount_cell2++;
            }
        }
    }

    free_1D<int>(sendbuf);
    free_1D<double>(sendcoor);
    free_1D<double>(recvcoor);
}
