#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
#include "pinit_data.h"
#include "pdistri_pro.h"
#include "pupdate_squares.h"


/*
 * Program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])
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
  // New variables required
  int rank, size;
  int left,right,up,down;
  MPI_Comm comm2d;

  /*
   *  Initialise MPI, compute the size
   */
  mp_start(&size);

  /*
   *  Cartesian topology
   */
  mp_cart(&comm2d, &M, &N, size, L);

  /*
   *  Find 4 neighbours of every squares for halo swaps
   */
  mp_find_neighbours(&rank, comm2d, &left, &right, &down, &up);

  /*
   *  Allocate memory dynamically
   */
  map = (int **) arralloc(sizeof(int), 2, L, L);
  smallmap  = (int**)arralloc(sizeof(int), 2, M, N);
  old = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);
  new = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);

  /*
   *  Initialize array map by master process
   */
  init_map(seed, rank, rho, map, L);

  /*
   *  Master process distributes array smallmap to others
   */
  mp_scatter_pro(map, smallmap, L, comm2d, M, N, rank);

  /*
   * Initialise the old array: copy the MxN array smallmap to the centre of
   * old, and set the halo values to zero.
   */
  init_old(smallmap, old, M, N);

  /*
   *  Update squares until there is no change between steps
   */
  update_squares(M, N, L, old, new, left, right, up, down, comm2d, rank);

  /*
   *  Copy the centre of old, excluding the halos, into smallmap
   */
  final_suqares(M, N, smallmap, old);

  /*
   *  Master process collects data back into array map
   */
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

  /*
   *  Finalize MPI
   */
  mp_stop();

  return 0;
}
