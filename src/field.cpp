#include "field.h"
#include "function.h"
#include "molecule.h"
#include "constants.h"

Field::Field(class Frame* iframe)
{
    fieldnum = NULL;
    fieldv = NULL;
    fieldT = NULL;

    fieldnum_all = NULL;
    fieldv_all = NULL;
    fieldT_all = NULL;

    frame = iframe;
}

Field::~Field()
{
    if(fieldnum)  free_2D<int>(fieldnum);
    if(fieldv)    free_2D<double>(fieldv);
    if(fieldT)    free_1D<double>(fieldT);

    if(fieldnum_all)  free_2D<int>(fieldnum);
    if(fieldv_all)    free_2D<double>(fieldv);
    if(fieldT_all)    free_1D<double>(fieldT);
}

void Field::cleanup_field()
{
    if(fieldnum)  free_2D<int>(fieldnum);
    if(fieldv)    free_2D<double>(fieldv);
    if(fieldT)    free_1D<double>(fieldT);

    if(fieldnum_all)  free_2D<int>(fieldnum_all);
    if(fieldv_all)    free_2D<double>(fieldv_all);
    if(fieldT_all)    free_1D<double>(fieldT_all);
}

void Field::cal_field(FILE* outfield)   
{
    int* lfield = frame->input->field_pars->lfield;
    if(!lfield[0] && !lfield[1] && !lfield[2]) return;

    int myid = frame->myid;
    cal_field_local();
    comm_field();
    if(myid==0) output_field(outfield);
    cleanup_field();
}

void Field::cal_field_local()
{
    int* lfield = frame->input->field_pars->lfield;
    if(!lfield[0] && !lfield[1] && !lfield[2]) return;

    int myid = frame->myid;
    int* fieldmesh_global = frame->input->field_pars->field_mesh;
    int nmolty = frame->input->sys_pars->nmolty;
    int* procgrid = frame->procgrid;
    double* boxlen = frame->domain->boxlen;
    int nlocal = frame->atom->nlocal;
    double** x = frame->atom->x;
    double** v = frame->atom->v;
    int* type = frame->atom->type;
    double* sublo = frame->domain->sublo;
    Molecule* all_molecule = frame->input->sys_pars->all_molecule;

    nmesh = 1;
    for(int i=0;i<3;i++)
    {
        fieldmesh[i] = fieldmesh_global[i]/procgrid[i];
        fieldlen[i] = boxlen[i]/procgrid[i];
        meshsize[i] = fieldlen[i]/fieldmesh[i];
        nmesh *= fieldmesh[i];
    }

    allocate_2D<int>(fieldnum,nmolty,nmesh,"fieldnum");
    for(int i=0;i<nmolty;i++)
    {
        for(int j=0;j<nmesh;j++)
        {
            fieldnum[i][j] = 0;
        }
    }
 
    double* fieldmass = NULL;
    int* fieldid = NULL;   
    if(lfield[1] || lfield[2])
    {
        allocate_2D<double>(fieldv,3,nmesh,"fieldv");
        allocate_1D<double>(fieldmass,nmesh,"fieldmass");

        for(int i=0;i<3;i++)
        {
            for(int j=0;j<nmesh;j++)
            {
                fieldv[i][j] = 0.;
            }
        }

        for(int i=0;i<nmesh;i++)
        {
            fieldmass[i] = 0.;
        }

        if(lfield[1])
        {
            allocate_1D<double>(fieldT,nmesh,"fieldT");
            allocate_1D<int>(fieldid,nlocal,"fieldid");  //remember
            for(int i=0;i<nmesh;i++)
            {
                fieldT[i] = 0.;
            }
        }
    }

    //calculate field properties
    for(int i=0;i<nlocal;i++)
    {
        int molty = type[i];
        int id_3D[3],id_1D;
        for(int j=0;j<3;j++)
        {
            id_3D[j] = static_cast<int>((x[i][j]-sublo[j])/meshsize[j]);
            id_pbc(id_3D[j],fieldmesh[j]);
        }
        id_1D = mapto1D(id_3D,fieldmesh);        
        fieldnum[molty-1][id_1D]++;
        if(fieldid) fieldid[i] = id_1D;
        if(fieldv)
        {
            double beadmass = all_molecule[molty-1].beadmass[0];
            for(int j=0;j<3;j++)
            {
                fieldv[j][id_1D] += beadmass * v[i][j]; 
            }
            fieldmass[id_1D] += beadmass;
        }
    }


    if(fieldv)
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<nmesh;j++)
            {
                double mass = fieldmass[j];
                if(mass > 0)
                {
                   fieldv[i][j] /= mass;
                }
            }
        }
    }

    if(fieldT)
    {
        for(int i=0;i<nlocal;i++)
        {
            int molty = type[i];
            int id_1D = fieldid[i];
            double vx_sub, vy_sub, vz_sub;
            double beadmass = all_molecule[molty-1].beadmass[0];
            vx_sub = v[i][0] - fieldv[0][id_1D];
            vy_sub = v[i][1] - fieldv[1][id_1D];
            vz_sub = v[i][2] - fieldv[2][id_1D];              
            fieldT[id_1D] += beadmass * (vx_sub*vx_sub + vy_sub*vy_sub + vz_sub*vz_sub);
        }

        for(int i=0;i<nmesh;i++)
        {
            int dof = 0;
            for(int j=0;j<nmolty;j++)
            {
                dof += fieldnum[j][i] * 3;
            }
            if(dof > 0)
            {
                double tfactor = constants::mvv2e /(dof*constants::boltz);
                fieldT[i] *= tfactor;
            }
        }
    } 

    if(fieldmass) free_1D<double>(fieldmass);
    if(fieldid)   free_1D<int>(fieldid);
}

void Field::comm_field()
{
    MPI_Datatype recvarray, resized_recvarray;
    int recvstarts[3],recvsubsizes[3],recvsizes[3];
    int myid = frame->myid;
    int* mycor = frame->mycor;
    int nproc = frame->nproc;
    int nmolty = frame->input->sys_pars->nmolty;
    int* fieldmesh_global = frame->input->field_pars->field_mesh;

    nmesh_all = 1;
    for(int i=0;i<3;i++)
    {
        nmesh_all *= fieldmesh_global[i];
    }


    if(myid!=0)
    {
        allocate_2D<int>(fieldnum_all,1,1,"fieldnum_all");        
        if(fieldv)
            allocate_2D<double>(fieldv_all,1,1,"fieldv_all");
        if(fieldT)
            allocate_1D<double>(fieldT_all,1, "fieldT_all");
    }
    else
    {
        allocate_2D<int>(fieldnum_all,nmolty,nmesh_all,"fieldnum_all");
        if(fieldv)
            allocate_2D<double>(fieldv_all,3,nmesh_all,"fieldv_all");
        if(fieldT)
            allocate_1D<double>(fieldT_all,nmesh_all,"fieldT_all");
    }

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
        int ndisp = fieldmesh_global[i]/frame->procgrid[i];
        mydisp_3D[i] = mycor[i]*ndisp;
    }

    mydisp_1D = mapto1D(mydisp_3D,fieldmesh_global);

    //communicate fieldnum
    MPI_Allgather(&mydisp_1D,1,MPI_INT,displs,1,MPI_INT,MPI_COMM_WORLD);

    for(int i=0;i<3;i++)
    {
        recvstarts[i] = 0;
        recvsubsizes[i] = fieldmesh[i];
        recvsizes[i] = fieldmesh_global[i];
    }

    MPI_Type_create_subarray(3,recvsizes,recvsubsizes,recvstarts,MPI_ORDER_C,MPI_INT,&recvarray);
    MPI_Type_commit(&recvarray);
    MPI_Type_create_resized(recvarray,0,sizeof(int),&resized_recvarray);
    MPI_Type_commit(&resized_recvarray);

    for(int i=0;i<nmolty;i++)
    {
        MPI_Gatherv(fieldnum[i],nmesh,MPI_INT,fieldnum_all[i],recvcounts,displs,resized_recvarray,0,MPI_COMM_WORLD);
    }
    MPI_Type_free(&recvarray);
    MPI_Type_free(&resized_recvarray);

    //communicate fieldv
    if(fieldv)
    {
        MPI_Type_create_subarray(3,recvsizes,recvsubsizes,recvstarts,MPI_ORDER_C,MPI_DOUBLE,&recvarray);
        MPI_Type_commit(&recvarray);
        MPI_Type_create_resized(recvarray,0,sizeof(double),&resized_recvarray);
        MPI_Type_commit(&resized_recvarray);

        for(int i=0;i<3;i++)
        {
            MPI_Gatherv(fieldv[i],nmesh,MPI_DOUBLE,fieldv_all[i],recvcounts,displs,resized_recvarray,0,MPI_COMM_WORLD);
        }
        MPI_Type_free(&recvarray);
        MPI_Type_free(&resized_recvarray);
    }

    //communicate fieldT
    if(fieldT)
    {
        MPI_Type_create_subarray(3,recvsizes,recvsubsizes,recvstarts,MPI_ORDER_C,MPI_DOUBLE,&recvarray);
        MPI_Type_commit(&recvarray);
        MPI_Type_create_resized(recvarray,0,sizeof(double),&resized_recvarray);
        MPI_Type_commit(&resized_recvarray);

        MPI_Gatherv(fieldT,nmesh,MPI_DOUBLE,fieldT_all,recvcounts,displs,resized_recvarray,0,MPI_COMM_WORLD);
        MPI_Type_free(&recvarray);
        MPI_Type_free(&resized_recvarray);
    }

    free_1D<int>(displs);
    free_1D<int>(recvcounts);
}

void Field::output_field(FILE* outfield)
{
    int nmolty = frame->input->sys_pars->nmolty;

    for(int i=0;i<nmesh_all;i++)
    {
        for(int j=0;j<nmolty;j++)
        {
            fprintf(outfield,"%15d ",fieldnum_all[j][i]);
        }
        if(fieldv_all) 
        {
            for(int j=0;j<3;j++)
            {
                fprintf(outfield,"%15.5e ", fieldv_all[j][i]*1e5); //convert to m/s
            }
        }
        if(fieldT_all)
        {
            fprintf(outfield,"%15.5e", fieldT_all[i]);
        }
        fprintf(outfield,"\n");
    }
}
