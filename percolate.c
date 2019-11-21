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
   *  Define the main arrays for the simulation
   */
  int L, M, N;
  int** old;
  int** new;


  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */
  int** map;

  /*
   *  Create a new M*N array called smallmap without any halos
   */
  int** smallmap;

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  //int step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc;

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
  //init_arrays(map, smallmap, old, new, L, M, N);

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

  int tag = 0;
  MPI_Status status;
  //MPI_Gather(&smallmap[0][0], M*N, MPI_INT, &map[0][0], M*N, MPI_INT, 0, comm);
  //MPI_Ssend(&smallmap[0][0], M*N, MPI_INT, 0, tag, comm2d);
  /*
  MPI_Datatype rectype;
  MPI_Type_vector(M,N,N,MPI_INT,&rectype);
  MPI_Type_commit(&rectype);*/
  int i,j;
  if(rank != 0){
     MPI_Ssend(&smallmap[0][0], M*N, MPI_INT, 0, tag, comm2d);
  }else{

     int coord[2];

     MPI_Cart_coords(comm2d, rank, 2, coord);

      int x = coord[0] * M;
      int y = coord[1] * N;

      for (i = 0; i < M; ++i) {
          for (j = 0; j < N ; ++j) {
             map[x+i][y+j] = smallmap[i][j];
          }
      }
      int source;
      for (source = 1; source < size; source++){

          /*MPI_Cart_coords(comm2d, source, 2, coord);

          int x = coord[0] * M;
          int y = coord[1] * N;

          MPI_Recv(&map[x][y], 1, rectype, source, tag, comm2d, &status);*/
          MPI_Recv(&smallmap[0][0], M*N, MPI_INT, source, tag, comm2d, &status);

          MPI_Cart_coords(comm2d, source, 2, coord);

          int x = coord[0] * M;
          int y = coord[1] * N;

          for (i = 0; i < M; ++i) {
              for (j = 0; j < N ; ++j) {
                  map[x+i][y+j] = smallmap[i][j];
              }
          }

      }
  }

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */
  if(rank == 0){
    perc = 0;

    for (itop=0; itop < L; itop++)
      {
        if (map[itop][L-1] > 0)
	 {
	   for (ibot=0; ibot < L; ibot++)
	      {
	        if (map[itop][L-1] == map[ibot][0])
		  {
		    perc = 1;
		  }
	     }
	 }
     }

   if (perc != 0)
     {
       printf("percolate: cluster DOES percolate\n");
     }
   else
     {
       printf("percolate: cluster DOES NOT percolate\n");
      }

  /*
   *  Write the map to the file "map.pgm", displaying only the very
   *  largest cluster (or multiple clusters if exactly the same size).
   *  If the last argument here was 2, it would display the largest 2
   *  clusters etc.
   */

    //percwrite("map.pgm", map, 8);
    percwritedynamic("map.pgm", map, L, 8);
  }

  MPI_Finalize();

  return 0;
}
