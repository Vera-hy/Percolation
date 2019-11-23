#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>
#include <math.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Cartesian topology
 */
void mp_cart(int *comm2d, int *m, int *n, int size, int l, int npro[]){

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

    npro[0] = dims[0];
    npro[1] = dims[1];
    *m = (int) floor(l / npro[0]);
    *n = (int) floor(l / npro[1]);

}

/*
 *  Find 4 neighbours of every squares for halo swaps
 *  and determines the rank of the calling process in
 *  the communicator.
 */
void mp_find_neighbours(int *rank, int comm2d, int *left, int *right,
        int *down, int *up, int *m, int *n, int npro[], int l){
    int directioni, directionj, disp, i, j;
    directioni = 0;       // shift along the first index
    directionj = 1;       // shift along the second index
    disp = 1;            // Shift by 1

    mpCommrank(comm2d, rank);

    for (j = 0; j <= npro[1] - 1; j++) {
        int coord1[2];
        int id1;
        coord1[0] = npro[0] - 1;
        coord1[1] = j;
        MPI_Cart_rank(comm2d, coord1, &id1);
        if (*rank == id1){
            *m = *m + (l - (npro[0] * (*m)));
        }
    }

    for (i = 0; i <= npro[0] - 1; i++) {
        int coord2[2];
        int id2;
        coord2[0] = i;
        coord2[1] = npro[1] - 1;
        MPI_Cart_rank(comm2d, coord2, &id2);
        if (*rank == id2){
            *n = *n + (l - (npro[1] * (*n)));
        }
    }

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


