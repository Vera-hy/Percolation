#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

/*
 * Copy the values of map to smallmap.
 */
void mp_scatter_pro(int** map, int** smallmap, int l,
                    int comm2d, int m, int n, int rank){
    int i, j;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n ; ++j) {
            smallmap[i][j] = map[i][j];
        }
    }
}
