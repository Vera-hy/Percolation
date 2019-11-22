#include <stdio.h>
#include <stdlib.h>

#include "percolate.h"

/*
 *  Set the value of size.
 */
void mp_start(int *size){
    *size = 1;
}

/*
 * @return 0 means jump out of this function successfully.
 */
int mp_stop(){
    return 0;
}

/*
 * Set the value of m, n and comm2d.
 */
void mp_cart(int *comm2d, int *m, int *n, int size, int l){
    *m = l;
    *n = l;
    *comm2d = 0;
}

/*
 * Set the value of rank, left, right, up, down.
 */
void mp_find_neighbours(int *rank, int comm2d, int *left,
                        int *right, int *down, int *up){
    *rank = 0;
    *left = 0;
    *right = 0;
    *up = 0;
    *down = 0;
}

/*
 *  Initialise map with density rho. Zero indicates rock, a positive
 *  value indicates a hole. For the algorithm to work, all the holes
 *  must be initialised with a unique integer.
 */
void init_map(int seed, int rank, double rho, int** map, int l){
    rinit(seed);
    int nhole, i ,j;
    double r;
    nhole = 0;

    for (i=0; i < l; i++)
    {
        for (j=0; j < l; j++)
        {
            r=uni();

            if(r < rho)
            {
                map[i][j] = 0;
            }
            else
            {
                nhole++;
                map[i][j] = nhole;
            }
        }
    }

    printf("percolate: rho = %f, actual density = %f\n",
           rho, 1.0 - ((double) nhole)/((double) l*l) );
}

