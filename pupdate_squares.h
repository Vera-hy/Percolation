#ifndef PUPDATE_SQUARES_H
#define PUPDATE_SQUARES_H

void init_old(int** smallmap, int** old, int M, int N);
void update_squares(int M, int N, int L, int** old, int** new, int left, int right,
                    int up, int down, MPI_Comm comm2d, int rank);
void final_suqares(int M, int N, int** smallmap, int** old);

#endif
