#ifndef PDISTRI_PRO_H
#define PDISTRI_PRO_H
void scatter_pro(int** map, int** smallmap, int L,
                 MPI_Comm comm2d, int M, int N, int rank);
#endif
