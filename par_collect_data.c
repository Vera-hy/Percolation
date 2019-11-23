#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>
#include <math.h>

#include "percolate.h"
#include "mplib.h"

/*
 *  Master process collects data back into array map.
 */
void mp_collect_data(int rank, int** smallmap, int m, int n, int comm2d,
        int size, int** map, int npro[], int l){
    int tag = 0;
    MPI_Status status;
    int i,j;

    if(rank != 0){
        //printf("rank = %d, m = %d, n = %d\n", rank, m, n);
        mpSsend(&smallmap[0][0], m*n, MPI_INT, 0, tag, comm2d);

    }else{

        for (i = 0; i < m; ++i) {
            for (j = 0; j < n; ++j) {
                map[i][j] = smallmap[i][j];
            }
        }

        int source;
        for (source = 1; source < size; source++){

            for (j = 0; j <= npro[1] - 1; j++) {
                int coord1[2];
                int id1;
                coord1[0] = npro[0] - 1;
                coord1[1] = j;
                MPI_Cart_rank(comm2d, coord1, &id1);
                if (source == id1){
                    m = m + (l - (npro[0] * m));
                }
            }

            for (i = 0; i <= npro[0] - 1; i++) {
                int coord2[2];
                int id2;
                coord2[0] = i;
                coord2[1] = npro[1] - 1;
                MPI_Cart_rank(comm2d, coord2, &id2);
                if (source == id2){
                    n = n + (l - (npro[1] * n));
                }
            }

            printf("source = %d, m = %d, n = %d\n", source, m, n);
           // printf("hi\n");
            mpRecv(&smallmap[0][0], m*n, MPI_INT, source, tag, comm2d, &status);

            m = (int) floor(l / npro[0]);
            n = (int) floor(l / npro[1]);

            int coord[2];
            mpCartcoords(comm2d, source, coord);

            int x = coord[0] * m;
            int y = coord[1] * n;

            for (j = 0; j <= npro[1] - 1; j++) {
                int coord5[2];
                int id5;
                coord5[0] = npro[0] - 1;
                coord5[1] = j;
                MPI_Cart_rank(comm2d, coord5, &id5);
                if (source == id5){
                    m = m + (l - (npro[0] * m));
                }
            }

            for (i = 0; i <= npro[0] - 1; i++) {
                int coord6[2];
                int id6;
                coord6[0] = i;
                coord6[1] = npro[1] - 1;
                MPI_Cart_rank(comm2d, coord6, &id6);
                if (source == id6){
                    n = n + (l - (npro[1] * n));
                }
            }

            for (i = 0; i < m; ++i) {
                for (j = 0; j < n ; ++j) {
                    map[x+i][y+j] = smallmap[i][j];
                }
            }
            m = (int) floor(l / npro[0]);
            n = (int) floor(l / npro[1]);

        }
    }
}
