#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "percolate.h"

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

    int step, oldval, newval, nchange, printfreq,i ,j;
    double tstart, tstop;

    printfreq = 100;
    step = 1;
    nchange = 1;

    tstart = MPI_Wtime();
    while (nchange > 0)
    {
        nchange = 0;

        for (j=1; j <= n; j++)
        {
            old[0][j]   = old[m][j];
            old[m+1][j] = old[1][j];
        }

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
        if (step % printfreq == 0)
        {
            printf("percolate: number of changes on step %d is %d\n",
                   step, nchange);
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
         * Calculate the average value of array map in each step
         * and print it to screen.
         */
        float map_average;
        int map_sum = 0;
        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                map_sum += old[i][j];
            }
        }
        map_average = map_sum / (l * l);
        if (step % printfreq == 0) {
            printf("percolate: the average of the map array on step %d "
                   "is %6.2f\n", step, map_average);
        }

        step++;
    }
    tstop = MPI_Wtime();
    printf("\nAverage time taken per step was  %f seconds\n", (tstop-tstart)/step);

}

/*
 * Copy the value of of array old back into the array smallmap.
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
