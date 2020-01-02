#include<cstdio>
#include<cstdlib>
#include<mpi.h>

#include "frame.h"

int main(int argc, char** argv)
{
    int myid;
    FILE* outbubble=NULL;
    FILE* outshell=NULL;
    FILE* outfield=NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(argc!=2)
    {
        printf("Incorrect number of arguments\n");
        exit(1);
    }

    Frame* frame = new Frame(argv[1]);
    frame->init();   //read input parameter, create subclass and cartisian communicator

    int* lbubble = frame->input->bubble_pars->lbubble;
    int* lshell = frame->input->shell_pars->lshell;
    int* lfield = frame->input->field_pars->lfield;

    if(myid==0)
    {
        if(lbubble[0] || lbubble[1]|| lbubble[2])
            outbubble = fopen("vdetail.dat","w");
        if(lshell[0] || lshell[1] || lshell[2])
            outshell = fopen("shell.dat","w");
        if(lfield[0] || lfield[1] || lfield[2])
            outfield = fopen("field.dat","w"); 
    }

    int maxframe = frame->input->file_pars->maxframe;
    double read_count = 0., bubble_count = 0., shell_count = 0., field_count = 0.;  //timer

    for(int iframe=0; iframe<maxframe; iframe++)
    {
        read_count -= MPI_Wtime();
        frame->read_head();
        frame->read_atom();
        printf("finish reading atom\n");
        read_count += MPI_Wtime();

        bubble_count -= MPI_Wtime();
        frame->bubble->cal_bubble(outbubble);
        printf("finish cal bubble\n");
        bubble_count += MPI_Wtime();

        shell_count -= MPI_Wtime();
        frame->shell->cal_shell(outshell);
        printf("finish cal shell\n");
        shell_count += MPI_Wtime();

        field_count -= MPI_Wtime();
        frame->field->cal_field(outfield);
        printf("finish cal field\n");
        field_count += MPI_Wtime();

        frame->atom->cleanup_atom();
        if(myid==0) printf("%d\n",iframe);
    }

    if(myid==0)
    {
        printf("Timing decomposition (min):\n");
        printf("Read: %10.3e\n",read_count/60.0);
        printf("Bubble: %10.3e\n",bubble_count/60.0);
        printf("Shell: %10.3e\n",shell_count/60.0);
        printf("Field: %10.3e\n",field_count/60.0);
    }

    delete frame;

    if(outbubble) fclose(outbubble);
    if(outshell)  fclose(outshell);
    if(outfield)  fclose(outfield);

    MPI_Finalize();
    return 0;
}
