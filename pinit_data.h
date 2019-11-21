#ifndef PINIT_DATA_H
#define PINIT_DATA_H

void mp_start(int *size);
void mp_stop();
void mp_cart(MPI_Comm *comm2d, int *M, int *N, int size, int L);
void mp_find_neighbours(int *rank, MPI_Comm comm2d, int *left,
                     int *right, int *down, int *up);
//void init_arrays(int** map, int** smallmap, int** old, int** new, int L,
//                 int M, int N);
void init_map(int seed, int rank, double rho, int** map, int L);

#endif
