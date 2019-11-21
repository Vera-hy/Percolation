#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

void mp_collect_data(int rank, int** smallmap, int M, int N, MPI_Comm comm2d,
        int size, int** map){
    int tag = 0;
    MPI_Status status;
    int i,j;

    if(rank != 0){
        MPI_Ssend(&smallmap[0][0], M*N, MPI_INT, 0, tag, comm2d);
    }else{

        int coord[2];

        MPI_Cart_coords(comm2d, rank, 2, coord);

        int x = coord[0] * M;
        int y = coord[1] * N;

        for (i = 0; i < M; ++i) {
            for (j = 0; j < N ; ++j) {
                map[x+i][y+j] = smallmap[i][j];
            }
        }
        int source;
        for (source = 1; source < size; source++){

            MPI_Recv(&smallmap[0][0], M*N, MPI_INT, source, tag, comm2d, &status);

            MPI_Cart_coords(comm2d, source, 2, coord);

            int x = coord[0] * M;
            int y = coord[1] * N;

            for (i = 0; i < M; ++i) {
                for (j = 0; j < N ; ++j) {
                    map[x+i][y+j] = smallmap[i][j];
                }
            }

        }
    }
}
