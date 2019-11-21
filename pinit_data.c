#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

/*
 *  Initialise MPI, compute the size
 */
void mp_start(int *size){

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, size);
}

/*
 *  Finalize MPI, compute the size
 */
void mp_stop(){
    MPI_Finalize();
}

/*
 *  Cartesian topology
 */
void mp_cart(MPI_Comm *comm2d, int *M, int *N, int size, int L){

    int reorder;
    int dims[ndims];
    int period[ndims];
    dims[0] = 0;
    dims[1] = 0;
    period[0] = TRUE;    // TRUE, Cyclic
    period[1] = FALSE;    // FALSE, Not Cyclic
    reorder = FALSE;      // FALSE

    MPI_Dims_create(size,ndims,dims);//dims on return 为每个维度P的数量
    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,period,reorder,comm2d);

    *M = L / dims[0];
    *N = L / dims[1];


//printf("My rank is: %d, my left is: %d, my right is:%d, my up is: %d, my down is: %d\n",
//         rank, left, right, up, down);
}

void mp_find_neighbours(int *rank, MPI_Comm comm2d, int *left,
                     int *right, int *down, int *up){
    int directioni, directionj, disp;
    directioni = 0;       // shift along the first index
    directionj = 1;       // shift along the second index
    disp = 1;            // Shift by 1

    MPI_Comm_rank(comm2d,rank);
    MPI_Cart_shift(comm2d,directioni,disp,left,right);
    MPI_Cart_shift(comm2d,directionj,disp,down,up);
}
/*
void init_arrays(int** map, int** smallmap, int** old, int** new, int L,
        int M, int N){
    map = (int **) arralloc(sizeof(int), 2, L, L);
    smallmap  = (int**)arralloc(sizeof(int), 2, M, N);//allocate memory
    old = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);//allocate memory
    new = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);//allocate memory
}*/


/*
 *  Initialise the generator and initialise map with density rho.
 *  Zero indicates rock, a positive value indicates a hole.
 *  For the algorithm to work, all the holes must be initialised with
 *   a unique integer. initialise map on the master process only
 */
void init_map(int seed, int rank, double rho, int** map, int L){

    int i, j, nhole;
    double r;

    rinit(seed);

    if(rank == 0){
        nhole = 0;

        for (i=0; i < L; i++)
        {
            for (j=0; j < L; j++)
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
               rho, 1.0 - ((double) nhole)/((double) L*L) );
    }
}


