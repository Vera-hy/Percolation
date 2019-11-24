#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"

#include "perclib.h"
#include "mplib.h"

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
  int l, m, n;
  int seed;
  double rho;

  if (argc != 3)
    {
      printf("Usage: percolate <seed> <l>\n");
      return 1;
    }

  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */
  rho = 0.411;

  /*
   *  Set the random number seed and length of the map array l
   */
  seed = atoi(argv[1]);

  l = atoi(argv[2]);

  // Define variables for parallel program
  int rank, size;
  int left,right,up,down;
  int comm2d;
  int npro[ndims];

  /*
   *  Initialise MPI, compute the size
   */
  mp_start(&size);

  /*
   *  Cartesian topology
   */
  mp_cart(&comm2d, &m, &n, size, l, npro);

  /*
   *  Find 4 neighbours of every squares for halo swaps
   *  and determines the rank of the calling process in
   *  the communicator.
   */
  //mp_find_neighbours(&rank, comm2d, &left, &right, &down, &up);
  mp_find_neighbours(&rank, comm2d, &left, &right, &down, &up,
            &m, &n, npro, l);
  //printf("rank = %d, m = %d, n = %d\n", rank, m,n);
  if (rank == 0){
      printf("npro0 = %d, npro1 = %d\n", npro[0],npro[1]);
      printf("percolate: params are l = %d, rho = %f, seed = %d\n", l, rho, seed);
  }

  /*
   *  Allocate memory dynamically
   */
  map = (int **) arralloc(sizeof(int), 2, l, l);
  if(rank != 0){
      smallmap  = (int**)arralloc(sizeof(int), 2, m, n);
      old = (int**)arralloc(sizeof(int), 2, m + 2, n + 2);
      new = (int**)arralloc(sizeof(int), 2, m + 2, n + 2);
  }
  if(rank == 0){
      int tempm = m + (l - (npro[0] * m));
      int tempn = n + (l - (npro[1] * n));
      smallmap  = (int**)arralloc(sizeof(int), 2, tempm, tempn);
      old = (int**)arralloc(sizeof(int), 2, tempm + 2, tempn + 2);
      new = (int**)arralloc(sizeof(int), 2, tempm + 2, tempn + 2);
  }


  /*
   *  Initialize array map by master process
   */
  init_map(seed, rank, rho, map, l);

  /*
   *  Master process distributes array map to others' array smallmap.
   */
  mp_scatter_pro(map, smallmap, l, comm2d, m, n, rank, npro);

  /*
   * Initialise the old array: copy the MxN array smallmap to the centre of
   * old, and set the halo values to zero.
   */
  init_old(smallmap, old, m, n);

  /*
   *  Update squares until there is no change between steps
   */
  update_squares(m, n, l, old, new, left, right, up, down, comm2d, rank, npro);

  /*
   *  Copy the centre of old, excluding the halos, back into array smallmap
   */
  final_suqares(m, n, smallmap, old);

  /*
   *  Master process collects data back into array map
   */
  mp_collect_data(rank, smallmap, m, n, comm2d, size, map, npro, l);


  if(rank == 0){
  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */
   test_perc(l, map);

  /*
   *  Write the map to the file "map.pgm", displaying only the very
   *  largest cluster (or multiple clusters if exactly the same size).
   *  If the last argument here was 8, it would display the largest 8
   *  clusters etc.
   */
    percwritedynamic("map.pgm", map, l, 8);
  }

  /*
   *  Finalize MPI
   */
  mp_stop();

  return 0;
}
