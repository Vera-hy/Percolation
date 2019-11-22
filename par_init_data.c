#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Cartesian topology
 */
void mp_cart(int *comm2d, int *m, int *n, int size, int l){

    int reorder;
    int dims[ndims];
    int period[ndims];
    dims[0] = 0;
    dims[1] = 0;
    period[0] = TRUE;    // TRUE, Cyclic
    period[1] = FALSE;    // FALSE, Not Cyclic
    reorder = FALSE;      // FALSE

    mpDimscreate(size, dims);//dims on return is the number of process in each dimension
    mpCartcreate(dims, period, reorder, comm2d);

    *m = l / dims[0];
    *n = l / dims[1];

}

/*
 *  Find 4 neighbours of every squares for halo swaps
 *  and determines the rank of the calling process in
 *  the communicator.
 */
void mp_find_neighbours(int *rank, int comm2d, int *left,
                     int *right, int *down, int *up){
    int directioni, directionj, disp;
    directioni = 0;       // shift along the first index
    directionj = 1;       // shift along the second index
    disp = 1;            // Shift by 1

    mpCommrank(comm2d, rank);

    mpCartshift(comm2d,directioni,disp,left,right);
    mpCartshift(comm2d,directionj,disp,down,up);
}

/*
 *  Initialise the generator and initialise map with density rho.
 *  Zero indicates rock, a positive value indicates a hole.
 *  For the algorithm to work, all the holes must be initialised with
 *   a unique integer. initialise map on the master process only.
 */
void init_map(int seed, int rank, double rho, int** map, int l){

    int i, j, nhole;
    double r;

    rinit(seed);

    if(rank == 0){
        nhole = 0;

        for (i=0; i < l; i++)
        {
            for (j=0; j < l; j++)
            {
                r=uni();

                if(r < rho)
                {
                    map[i][j] = 0;
                }
                else
                {
                    nhole++;
                    map[i][j] = nhole;
                }
            }
        }

        printf("percolate: rho = %f, actual density = %f\n",
               rho, 1.0 - ((double) nhole)/((double) l*l) );
    }
}


