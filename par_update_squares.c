#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "mplib.h"

void calc_max(int i, int j, int *nchange, int** old, int** new);

void overlap_swap_calc(int m, int n, int** old, int right, int left, int up,
        int down, int comm2d, int *nchange, int** new, int l, int npro[], int rank);

int calc_nchange(int nchange, int *all_nchange, int comm2d);

int calc_avemap(int m, int n, int** old, int comm2d, int l);

/*
 * Initialise the old array: copy the mxn array smallmap to the centre of
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
        int up, int down, int comm2d, int rank, int npro[]){

    int step, nchange, printfreq, all_nchange, i ,j;
    double tstart, tstop;
    //int maxstep = 40000;
    printfreq = 100;

    step = 1;
    nchange = 1;
    all_nchange = 1;

    tstart = gettime();
    //while (step <= maxstep)
    while (all_nchange > 0)
    {
        nchange = 0;
        all_nchange = 0;

        /*
         * Doing calculation that does not require the
         * communicated halo data at the same time as
         * the halo is being sent using non-blocking routines.
         */
        overlap_swap_calc(m, n, old, right, left, up,
                down, comm2d, &nchange, new, l, npro, rank);

        /*
         * Calculations involving the halos are done
         * separately after the halos have arrived.
         */
        for (i = 1; i <= m; i++) {
            calc_max(i, 1, &nchange, old, new);
            calc_max(i, n, &nchange, old, new);
        }
        for (j = 2; j <= n-1; j++) {
            calc_max(1, j, &nchange, old, new);
            calc_max(m, j, &nchange, old, new);
        }

        /*
         * Calculate and print out the number of changes of all the processes
         * in each step.
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

        /*
         * Calculate and print out the average value of array map in each step.
         */
        float map_average;
        map_average = calc_avemap(m, n, old, comm2d, l);

        if( rank == 0) {
            if (step % printfreq == 0) {
                printf("percolate: the average of the map array on step %d "
                       "is %6.2f\n", step, map_average);
            }
        }

        step++;
    }

    tstop = gettime();

    if(rank == 0){
        printf("\npercolate: parallel time taken in %d steps was %f seconds\n", step, tstop-tstart);
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
void final_squares(int m, int n, int** smallmap, int** old){

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
 * Set new[i][j] to be the maximum value of old[i][j]
 * and its four nearest neighbours.
 */
void calc_max(int i, int j, int *nchange, int** old, int** new){

    int oldval, newval;
    oldval = old[i][j];
    newval = oldval;

    if (oldval != 0)
    {
        if (old[i-1][j] > newval) newval = old[i-1][j];
        if (old[i+1][j] > newval) newval = old[i+1][j];
        if (old[i][j-1] > newval) newval = old[i][j-1];
        if (old[i][j+1] > newval) newval = old[i][j+1];

        if (newval != oldval)
        {
            ++ *nchange;
        }
    }

    new[i][j] = newval;

}

/*
* Doing calculation that does not require the communicated halo data at
* the same time as the halo is being sent using non-blocking routines.
*/
void overlap_swap_calc(int m, int n, int** old, int right, int left,
        int up, int down, int comm2d, int *nchange, int** new, int l, int npro[], int rank){

    MPI_Status recv_status[4], send_status[4];
    MPI_Request send_requests[4], recv_requests[4];
    int tag[4] = {1, 2, 3, 4};

    MPI_Datatype halo_rowtype;
    if(rank == 0){
        int tempn = n + (l - (npro[1] * n));
        mpVector(m, 1, tempn+2, MPI_INT, &halo_rowtype);
        mpTypecommit(&halo_rowtype);
    } else{
        mpVector(m, 1, n+2, MPI_INT, &halo_rowtype);
        mpTypecommit(&halo_rowtype);
    }

    mpIssend(&old[m][1], n, MPI_INT, right, tag[0], comm2d, &send_requests[0]);
    mpIssend(&old[1][1], n, MPI_INT, left, tag[1], comm2d, &send_requests[1]);
    mpIssend(&old[1][n], 1, halo_rowtype, up, tag[2], comm2d, &send_requests[2]);
    mpIssend(&old[1][1], 1, halo_rowtype, down, tag[3], comm2d, &send_requests[3]);

    mpIrecv(&old[0][1], n, MPI_INT, left, tag[0], comm2d, &recv_requests[0]);
    mpIrecv(&old[m+1][1], n, MPI_INT, right, tag[1], comm2d, &recv_requests[1]);
    mpIrecv(&old[1][0], 1, halo_rowtype, down, tag[2], comm2d, &recv_requests[2]);
    mpIrecv(&old[1][n+1], 1, halo_rowtype, up, tag[3], comm2d, &recv_requests[3]);

    int i,j;
    for (i = 2; i < m; i++) {
        for (j = 2; j < n; j++) {
            calc_max(i, j, nchange, old, new);
        }
    }

    mpWaitall(4, recv_requests, recv_status);
    mpWaitall(4, send_requests, send_status);

}

/*
 * Calculate the number of changes of all the processes in each step
 * so that the program can stop when there is no change between steps.
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


