#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>
#include <math.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Master process distributes array map to others' array smallmap.
 */
void mp_scatter_pro(int** map, int** smallmap, int l,
        int comm2d, int m, int n, int rank, int npro[]){

    mpBcast(&map[0][0], l*l, comm2d);

    int coord[2], i, j, tempm, tempn;

    tempm = (int) floor(l / npro[0]);
    tempn = (int) floor(l / npro[1]);

    mpCartcoords(comm2d, rank, coord);
    int x = coord[0] * tempm;
    int y = coord[1] * tempn;

    for (i = 0; i < m; ++i) {
        for (j = 0; j < n ; ++j) {
            smallmap[i][j] = map[x+i][y+j];
        }
    }
    //printf("rank = %d, m = %d, n = %d\n", rank, m,n);
}

