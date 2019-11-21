#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

void swap_halos(int M, int N, int** old, int right, int left, int up, int down, MPI_Comm comm2d);
int calc_nchange(int nchange, int *anchange, MPI_Comm comm2d);
int calc_avemap(int M, int N, int** old, MPI_Comm comm2d, int L);

/*
* Initialise the old array: copy the MxN array smallmap to the centre of
* old, and set the halo values to zero.
*/
void init_old(int** smallmap, int** old, int M, int N){

    int i,j;
    for (i=1; i <= M; i++)
    {
        for (j=1; j <= N; j++)
        {
            old[i][j] = smallmap[i-1][j-1];
        }
    }

    for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
        old[i][0]   = 0;
        old[i][N+1] = 0;
    }

    for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
        old[0][j]   = 0;
        old[M+1][j] = 0;
    }
}


/*
 *  Update for a fixed number of iterations
 */
void update_squares(int M, int N, int L, int** old, int** new, int left, int right,
        int up, int down, MPI_Comm comm2d, int rank){

    int step, oldval, newval, nchange, printfreq, anchange,i ,j;
    //maxstep = 16*L;
    printfreq = 100;

    step = 1;
    nchange = 1;
    anchange = 1;

//while (step <= maxstep)
    while (anchange > 0)
    {
        nchange = 0;
        anchange = 0;
        swap_halos(M, N, old, right, left, up, down, comm2d);
        for (i=1; i<=M; i++)
        {
            for (j=1; j<=N; j++)
            {
                oldval = old[i][j];
                newval = oldval;
        /*
         * Set new[i][j] to be the maximum value of old[i][j]
         * and its four nearest neighbours
         */
                if (oldval != 0)
                {
                    if (old[i-1][j] > newval) newval = old[i-1][j];
                    if (old[i+1][j] > newval) newval = old[i+1][j];
                    if (old[i][j-1] > newval) newval = old[i][j-1];
                    if (old[i][j+1] > newval) newval = old[i][j+1];

                    if (newval != oldval)
                    {
                        ++nchange;
                    }
                }

                new[i][j] = newval;
            }
        }

        /*
         *  Report progress every now and then
         */

        if (step % printfreq == 0)
        {
            printf("percolate: number of changes on step %d is %d\n",
                   step, nchange);
        }
        anchange = calc_nchange(nchange, &anchange, comm2d);

        // printf("anchange = %d\n ",anchange);

        /*
         *  Copy back in preparation for next step, omitting halos
         */

        for (i=1; i<=M; i++)
        {
            for (j=1; j<=N; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        float map_average;
        map_average = calc_avemap(M, N, old, comm2d, L);

        if( rank == 0) {
            if (step % printfreq == 0) {
                printf("percolate: the average of the map array on step %d "
                       "is %f\n", step, map_average);
            }
        }

        step++;
    }

        /*
         *  We set a maximum number of steps to ensure the algorithm always
         *  terminates. However, if we hit this limit before the algorithm
         *  has finished then there must have been a problem (e.g. maxstep
         *  is too small)
         */

    /*if (nchange != 0)
    {
        printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
               maxstep);
    }*/
}
/*
 *  Copy the centre of old, excluding the halos, into smallmap
 */
void final_suqares(int M, int N, int** smallmap, int** old){

    int i, j;
    for (i=1; i<=M; i++)
    {
        for (j=1; j<=N; j++)
        {
            smallmap[i-1][j-1] = old[i][j];
        }
    }
}

void swap_halos(int M, int N, int** old, int right, int left,
        int up, int down, MPI_Comm comm2d){

    MPI_Status status1,status2,status3,status4;
    MPI_Request request1,request2,request3,request4;
    int tag1 = 1;
    int tag2 = 2;
    int tag3 = 3;
    int tag4 = 4;

    MPI_Issend(&old[M][1],N,MPI_INT,right,tag1,comm2d,&request1);
    MPI_Recv(&old[0][1],N,MPI_INT,left,tag1,comm2d,&status1);
    MPI_Wait(&request1,&status1);

    MPI_Issend(&old[1][1],N,MPI_INT,left,tag2,comm2d,&request2);
    MPI_Recv(&old[M+1][1],N,MPI_INT,right,tag2,comm2d,&status2);
    MPI_Wait(&request2,&status2);

    MPI_Datatype hrowtype;
    MPI_Type_vector(M,1,N+2,MPI_INT,&hrowtype);
    MPI_Type_commit(&hrowtype);

    MPI_Issend(&old[1][N],1,hrowtype,up,tag3,comm2d,&request3);
    MPI_Recv(&old[1][0],1,hrowtype,down,tag3,comm2d,&status3);
    MPI_Wait(&request3,&status3);

    MPI_Issend(&old[1][1],1,hrowtype,down,tag4,comm2d,&request4);
    MPI_Recv(&old[1][N+1],1,hrowtype,up,tag4,comm2d,&status4);
    MPI_Wait(&request4,&status4);
}

int calc_nchange(int nchange, int *anchange, MPI_Comm comm2d){

    MPI_Allreduce(&nchange,anchange,1,MPI_INT,MPI_SUM,comm2d);
    return *anchange;
}

int calc_avemap(int M, int N, int** old, MPI_Comm comm2d, int L){
    int smallmap_sum = 0;
    int map_sum;
    float map_average;
    int i, j;
    for (i=1; i<=M; i++)
    {
        for (j=1; j<=N; j++)
        {
            smallmap_sum += old[i][j];
        }
    }
    MPI_Allreduce(&smallmap_sum,&map_sum,1,MPI_INT,MPI_SUM,comm2d);
    map_average = map_sum / (L * L);
    return  map_average;
}


