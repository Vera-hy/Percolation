#ifndef SEQLIB_H
#define SEQLIB_H

/*
 *  Prototypes for functions of sequential version.
 */

/*
 *  Set the value of size.
 */
void mp_start(int *size);

/*
 * @return 0 means jump out of this function successfully.
 */
int mp_stop();

/*
 * Set the value of m, n and comm2d.
 */
void mp_cart(int *comm2d, int *m, int *n, int size, int l);

/*
 * Set the value of rank, left, right, up, down.
 */
void mp_find_neighbours(int *rank, int comm2d, int *left,
                        int *right, int *down, int *up);

/*
 *  Initialise map with density rho.
 */
void init_map(int seed, int rank, double rho, int** map, int l);

/*
 * Copy the values of map to smallmap.
 */
void mp_scatter_pro(int** map, int** smallmap, int l,
                    int comm2d, int m, int n, int rank);

/*
 * Initialise the old array: copy the MxN array smallmap to the centre of
 * old, and set the halo values to zero.
 */
void init_old(int** smallmap, int** old, int m, int n);

/*
 *  Update squares until there is no change between steps.
 */
void update_squares(int m, int n, int l, int** old, int** new, int left, int right,
int up, int down, int comm2d, int rank);

/*
 * Copy the value of of array old back into the array smallmap.
 */
void final_suqares(int m, int n, int** smallmap, int** old);

/*
 * Copy the values of smallmap to map.
 */
void mp_collect_data(int rank, int** smallmap, int m, int n, int comm2d,
                     int size, int** map);

#endif
