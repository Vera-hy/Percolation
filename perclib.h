#ifndef PARLIB_H
#define PARLIB_H

/*
 *  Prototypes for functions of parallel version and sequential version.
 */

/*
 *  Cartesian topology
 */
void mp_cart(int *comm2d, int size, int l, int npro[]);

/*
 *  Find 4 neighbours of every squares for halo swaps,
 *  determines the rank of the calling process in
 *  the communicator and decide the value of m and n
 *  for array smallmap.
 */
void mp_find_neighbours(int *rank, int comm2d, int *left, int *right,
                        int *down, int *up, int *m, int *n, int npro[], int l);

/*
 *  Initialize array map by master process.
 */
void init_map(int seed, int rank, double rho, int** map, int l);

/*
 *  Master process distributes array map to others' array smallmap.
 */
void mp_scatter_pro(int** map, int** smallmap, int l,
                    int comm2d, int m, int n, int rank, int npro[]);

/*
 * Initialise the old array: copy the MxN array smallmap to the centre of
 * old, and set the halo values to zero.
 */
void init_old(int** smallmap, int** old, int m, int n);

/*
 *  Update squares until there is no change between steps.
 */
void update_squares(int m, int n, int l, int** old, int** new, int left, int right,
int up, int down, int comm2d, int rank, int npro[]);

/*
 *  Copy the centre of old, excluding the halos, back into smallmap.
 */
void final_suqares(int m, int n, int** smallmap, int** old);

/*
 *  Master process collects data back into array map.
 */
void mp_collect_data(int rank, int** smallmap, int m, int n, int comm2d,
                     int size, int** map, int npro[], int l);

#endif
