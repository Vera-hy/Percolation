#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Master process distributes array map to others' array smallmap.
 */
void mp_scatter_pro(int** map, int** smallmap, int l,
        int comm2d, int m, int n, int rank){

    mpBcast(&map[0][0], l*l, comm2d);

    int coord[2];

    mpCartcoords(comm2d, rank, coord);

    int x = coord[0] * m;
    int y = coord[1] * n;
    int i,j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n ; ++j) {
            smallmap[i][j] = map[x+i][y+j];
        }
    }
}

