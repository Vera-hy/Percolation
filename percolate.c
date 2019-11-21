#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "pinit_data.h"
#include "pdistri_pro.h"
#include "pupdate_squares.h"


/*
 * Serial program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])//test for Clion
{
  /*
   *  Define the arrays for the simulation
   */
  int** old;
  int** new;
  int** map;
  int** smallmap;

  /*
   *  Variables that define the simulation
   */
  int L, M, N;
  int seed;
  double rho;

  if (argc != 3)
    {
      printf("Usage: percolate <seed> <L>\n");
      return 1;
    }

  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */

  rho = 0.411;

  /*
   *  Set the randum number seed and length of the map array L
   */

  seed = atoi(argv[1]);

  L = atoi(argv[2]);

  printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);


  /*
   *  initialise MPI, compute the size and rank, and check that size = P
   */

  int rank, size;
  int left,right,up,down;
  MPI_Comm comm2d;

  mp_start(&size);
  mp_cart(&comm2d, &M, &N, size, L);
  mp_find_neighbours(&rank, comm2d, &left, &right, &down, &up);

  /*
   *  Allocate memory
   */
  map = (int **) arralloc(sizeof(int), 2, L, L);
  smallmap  = (int**)arralloc(sizeof(int), 2, M, N);
  old = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);
  new = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);

  init_map(seed, rank, rho, map, L);

  mp_scatter_pro(map, smallmap, L, comm2d, M, N, rank);

  init_old(smallmap, old, M, N);
  update_squares(M, N, L, old, new, left, right, up, down, comm2d, rank);
  final_suqares(M, N, smallmap, old);

  mp_collect_data(rank, smallmap, M, N, comm2d, size, map);

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */
  if(rank == 0){
   test_perc(L, map);

  /*
   *  Write the map to the file "map.pgm", displaying only the very
   *  largest cluster (or multiple clusters if exactly the same size).
   *  If the last argument here was 2, it would display the largest 2
   *  clusters etc.
   */

    percwritedynamic("map.pgm", map, L, 8);
  }

  MPI_Finalize();

  return 0;
}
