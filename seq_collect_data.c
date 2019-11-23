#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

/*
 * Copy the values of smallmap to map.
 */
void mp_collect_data(int rank, int** smallmap, int m, int n, int comm2d,
                     int size, int** map, int npro[], int l){
    int i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n ; ++j) {
            map[i][j] = smallmap[i][j];
        }
    }
}
