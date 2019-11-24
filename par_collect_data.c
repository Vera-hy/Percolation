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

        mpSsend(&smallmap[0][0], m*n, MPI_INT, 0, tag, comm2d);

    }else{

        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                map[i][j] = smallmap[i][j];
            }
        }

        int source;
        int stride = n + (l - (npro[1] * n));

        for (source = 1; source < size; source++){
            for (j = 0; j <= npro[1] - 1; j++) {
                int coord1[2];
                int id1;
                coord1[0] = npro[0] - 1;
                coord1[1] = j;
                mpCartrank(comm2d, coord1, &id1);
                if (source == id1){
                    m = m + (l - (npro[0] * m));
                }
            }

            for (i = 0; i <= npro[0] - 1; i++) {
                int coord2[2];
                int id2;
                coord2[0] = i;
                coord2[1] = npro[1] - 1;
               mpCartrank(comm2d, coord2, &id2);
                if (source == id2){
                    n = n + (l - (npro[1] * n));
                }
            }

            MPI_Datatype vector_mn;
            mpVector(m, n, stride, MPI_INT, &vector_mn);
            mpTypecommit(&vector_mn);

            mpRecv(&smallmap[0][0], 1, vector_mn, source, tag, comm2d, &status);

            int tempm = (int) floor(l / npro[0]);
            int tempn = (int) floor(l / npro[1]);

            int coord[2];
            mpCartcoords(comm2d, source, coord);

            int x = coord[0] * tempm;
            int y = coord[1] * tempn;

            for (i = 0; i < m; i++) {
                for (j = 0; j < n ; j++) {
                    map[x+i][y+j] = smallmap[i][j];
                }
            }
            m = tempm;
            n = tempn;

        }
    }
}
