#ifndef PCOLLECT_DATA_H
#define PCOLLECT_DATA_H

void mp_collect_data(int rank, int** smallmap, int M, int N, MPI_Comm comm2d,
                  int size, int** map);

#endif
