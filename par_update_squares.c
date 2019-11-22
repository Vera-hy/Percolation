#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

void swap_halos(int m, int n, int** old, int right, int left, int up, int down, int comm2d);
int calc_nchange(int nchange, int *all_nchange, int comm2d);
int calc_avemap(int m, int n, int** old, int comm2d, int l);

/*
 * Initialise the old array: copy the MxN array smallmap to the centre of
 * old, and set the halo values to zero.
 */
void init_old(int** smallmap, int** old, int m, int n){

    int i,j;
    for (i=1; i <= m; i++)
    {
        for (j=1; j <= n; j++)
        {
            old[i][j] = smallmap[i-1][j-1];
        }
    }

    for (i=0; i <= m+1; i++)  // zero the bottom and top halos
    {
        old[i][0]   = 0;
        old[i][n+1] = 0;
    }

    for (j=0; j <= n+1; j++)  // zero the left and right halos
    {
        old[0][j]   = 0;
        old[m+1][j] = 0;
    }
}

/*
 *  Update squares until there is no change between steps.
 */
void update_squares(int m, int n, int l, int** old, int** new, int left, int right,
        int up, int down, int comm2d, int rank){

    int step, oldval, newval, nchange, printfreq, all_nchange,i ,j;
    //maxstep = 16*l;
    printfreq = 100;

    step = 1;
    nchange = 1;
    all_nchange = 1;

    //while (step <= maxstep)
    while (all_nchange > 0)
    {
        nchange = 0;
        all_nchange = 0;
        swap_halos(m, n, old, right, left, up, down, comm2d);
        for (i=1; i<=m; i++)
        {
            for (j=1; j<=n; j++)
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
        all_nchange = calc_nchange(nchange, &all_nchange, comm2d);

        if (step % printfreq == 0) {
            if (rank == 0) {
                printf("percolate: number of all changes on step %d is %d\n",
                       step, all_nchange);
            }
        }
        /*
         *  Copy back in preparation for next step, omitting halos
         */
        for (i=1; i<=m; i++)
        {
            for (j=1; j<=n; j++)
            {
                old[i][j] = new[i][j];
            }
        }

        float map_average;
        map_average = calc_avemap(m, n, old, comm2d, l);

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
 *  Copy the centre of old, excluding the halos, back into smallmap.
 */
void final_suqares(int m, int n, int** smallmap, int** old){

    int i, j;
    for (i=1; i<=m; i++)
    {
        for (j=1; j<=n; j++)
        {
            smallmap[i-1][j-1] = old[i][j];
        }
    }
}

/*
 * Every process swaps halos to get correct values of halos.
 */
void swap_halos(int m, int n, int** old, int right, int left,
        int up, int down, int comm2d){

    MPI_Status status1,status2,status3,status4;
    MPI_Request request1,request2,request3,request4;
    int tag1 = 1;
    int tag2 = 2;
    int tag3 = 3;
    int tag4 = 4;

    mpIssend(&old[m][1], n, MPI_INT, right, tag1, comm2d, &request1);
    mpRecv(&old[0][1], n, MPI_INT, left, tag1, comm2d, &status1);
    mpWait(&request1, &status1);

    mpIssend(&old[1][1], n, MPI_INT, left, tag2, comm2d, &request2);
    mpRecv(&old[m+1][1], n, MPI_INT, right, tag2, comm2d, &status2);
    mpWait(&request2, &status2);

    MPI_Datatype halo_rowtype;
    mpVector(m, 1, n+2, MPI_INT, &halo_rowtype);
    mpTypecommit(&halo_rowtype);

    mpIssend(&old[1][n], 1, halo_rowtype, up, tag3, comm2d, &request3);
    mpRecv(&old[1][0], 1, halo_rowtype, down, tag3, comm2d, &status3);
    mpWait(&request3, &status3);

    mpIssend(&old[1][1], 1, halo_rowtype, down, tag4, comm2d, &request4);
    mpRecv(&old[1][n+1], 1, halo_rowtype, up, tag4, comm2d, &status4);
    mpWait(&request4, &status4);
}

/*
 * Calculate the number of changes of all the processes so that the program
 * can stop when there is no change between steps.
 *
 * @return the number of changes of all the processes
 */
int calc_nchange(int nchange, int *all_nchange, int comm2d){

    mpgsum(&nchange, all_nchange, 1, comm2d);
    return *all_nchange;
}

/*
 * Calculate the average value of array map in each step.
 *
 * @return the average value of array map in each step
 */
int calc_avemap(int m, int n, int** old, int comm2d, int l){
    int smallmap_sum = 0;
    int map_sum;
    float map_average;
    int i, j;
    for (i=1; i<=m; i++)
    {
        for (j=1; j<=n; j++)
        {
            smallmap_sum += old[i][j];
        }
    }
    mpgsum(&smallmap_sum, &map_sum, 1, comm2d);
    map_average = map_sum / (l * l);
    return  map_average;
}


