#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

/*
 * Scatter the map from the master to all other processes
 */
void scatter_pro(int** map, int** smallmap, int L,
        MPI_Comm comm2d, int M, int N, int rank){

    MPI_Bcast (&map[0][0], L*L, MPI_INT, 0, comm2d);

    int coord[2];

    MPI_Cart_coords(comm2d, rank, 2, coord);

    int x = coord[0] * M;
    int y = coord[1] * N;
    int i,j;
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N ; ++j) {
            smallmap[i][j] = map[x+i][y+j];
        }
    }
}

