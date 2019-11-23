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

    for (j = 0; j <= npro[1] - 1; j++) {
        int coord1[2];
        int id1;
        coord1[0] = npro[0] - 1;
        coord1[1] = j;
        MPI_Cart_rank(comm2d, coord1, &id1);
        if (rank == id1){
            tempm = m;
            m = (int) floor(l / npro[0]);
        }
    }

    for (i = 0; i <= npro[0] - 1; i++) {
        int coord2[2];
        int id2;
        coord2[0] = i;
        coord2[1] = npro[1] - 1;
        MPI_Cart_rank(comm2d, coord2, &id2);
        if (rank == id2){
            tempn = n;
            n = (int) floor(l / npro[1]);
        }
    }

    mpCartcoords(comm2d, rank, coord);
    int x = coord[0] * m;
    int y = coord[1] * n;
    //int i,j;
    for (j = 0; j <= npro[1] - 1; j++) {
        int coord3[2];
        int id3;
        coord3[0] = npro[0] - 1;
        coord3[1] = j;
        MPI_Cart_rank(comm2d, coord3, &id3);
        if (rank == id3){
            m = tempm;
        }
    }

    for (i = 0; i <= npro[0] - 1; i++) {
        int coord4[2];
        int id4;
        coord4[0] = i;
        coord4[1] = npro[1] - 1;
        MPI_Cart_rank(comm2d, coord4, &id4);
        if (rank == id4){
            n = tempn;
        }
    }
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n ; ++j) {
            smallmap[i][j] = map[x+i][y+j];
        }
    }
    //printf("rank = %d, m = %d, n = %d\n", rank, m,n);
}

