#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Master process collects data back into array map.
 */
void mp_collect_data(int rank, int** smallmap, int m, int n, int comm2d,
        int size, int** map){
    int tag = 0;
    MPI_Status status;
    int i,j;

    if(rank != 0){
        mpSsend(&smallmap[0][0], m*n, MPI_INT, 0, tag, comm2d);
    }else{

        int coord[2];

        mpCartcoords(comm2d, rank, coord);

        int x = coord[0] * m;
        int y = coord[1] * n;

        for (i = 0; i < m; ++i) {
            for (j = 0; j < n ; ++j) {
                map[x+i][y+j] = smallmap[i][j];
            }
        }
        int source;
        for (source = 1; source < size; source++){

            mpRecv(&smallmap[0][0], m*n, MPI_INT, source, tag, comm2d, &status);

            mpCartcoords(comm2d, source, coord);

            int x = coord[0] * m;
            int y = coord[1] * n;

            for (i = 0; i < m; ++i) {
                for (j = 0; j < n ; ++j) {
                    map[x+i][y+j] = smallmap[i][j];
                }
            }

        }
    }
}
