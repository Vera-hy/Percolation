#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdarg.h>

#include "percolate.h"
//#include "arralloc.h"


/*
 * Serial program to test for percolation of a cluster.
 */

int main(int argc, char *argv[])//test for Clion
{
  /*
   *  Define the main arrays for the simulation
   */
  int L, M, N;
  int **map;
  int** old;
  int** new;
  int** smallmap;
  //int old[M+2][N+2], new[M+2][N+2];//动态分配？

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  //int map[L][L];

  /*
   *  Create a new M*N array called smallmap without any halos
   */

  //int smallmap[M][N];//动态分配？

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc;
  double r;

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
   *  Set the randum number seed and initialise the generator
   */

  seed = atoi(argv[1]);

  L = atoi(argv[2]);

  printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);

  map = (int **) arralloc(sizeof(int), 2, L, L);

  rinit(seed);

  /*
   *  initialise MPI, compute the size and rank, and check that size = P
   */

  MPI_Comm comm;
  MPI_Status status1,status2,status3,status4;
  MPI_Request request1,request2,request3,request4;
    
  comm = MPI_COMM_WORLD;
    
  int rank,size,tag1,tag2,tag3,tag4;

  tag1 = 1;
  tag2 = 2;
  tag3 = 3;
  tag4 = 4;

  // New variables required
  int comm2d,directioni,directionj,disp;
  int dims[ndims];
  int period[ndims];
  int reorder,left,right,up,down;

  MPI_Init(NULL, NULL);

  MPI_Comm_size(comm, &size);
  //MPI_Comm_rank(comm, &rank);

  /*
   *  Cartesian topology
   */ 
  dims[0] = 0;
  dims[1] = 0;
  period[0] = TRUE;    // TRUE, Cyclic
  period[1] = FALSE;    // FALSE, Not Cyclic
  reorder = FALSE;      // FALSE
  directioni = 0;       // shift along the first index
  directionj = 1;       // shift along the second index
  disp = 1;            // Shift by 1

  MPI_Dims_create(size,ndims,dims);//dims on return 为每个维度P的数量
  MPI_Cart_create(comm,ndims,dims,period,reorder,&comm2d);

  M = L / dims[0];
  N = L / dims[1];

  MPI_Comm_rank(comm2d,&rank);
  MPI_Cart_shift(comm2d,directioni,disp,&left,&right);
  MPI_Cart_shift(comm2d,directionj,disp,&down,&up);
  //printf("My rank is: %d, my left is: %d, my right is:%d, my up is: %d, my down is: %d\n",
   //         rank, left, right, up, down);

  smallmap  = (int**)arralloc(sizeof(int), 2, M, N);//allocate memory
  old = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);//allocate memory
  new = (int**)arralloc(sizeof(int), 2, M + 2, N + 2);//allocate memory

  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer. initialise map on the
   *  master process only
   */
  if(rank == 0){
    nhole = 0;

    for (i=0; i < L; i++)
     {
        for (j=0; j < L; j++)
	 {
	   r=uni();
	  
	   if(r < rho)
	     {
	       map[i][j] = 0;
	     }
	   else
	      {
	        nhole++;
	        map[i][j] = nhole;
	      }
	 }
      }

    printf("percolate: rho = %f, actual density = %f\n",
	   rho, 1.0 - ((double) nhole)/((double) L*L) );
  }

  /*
   * Scatter the map from the master to all other processes using
   * MPI_Scatter with sendbuf = map and recvbuf = smallmap.
   */

  //MPI_Scatter(&map[0][0], M*N, MPI_INT, &smallmap[0][0], M*N, MPI_INT, 0, comm);
    MPI_Bcast (&map[0][0], L*L, MPI_INT, 0, comm2d);

    int coord[2];

    MPI_Cart_coords(comm2d, rank, 2, coord);

    int x = coord[0] * M;
    int y = coord[1] * N;

    for (i = 0; i < M; ++i) {
        for (j = 0; j < N ; ++j) {
            smallmap[i][j] = map[x+i][y+j];
        }
    }
  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */

   for (i=1; i <= M; i++)
    {
      for (j=1; j <= N; j++)
	    {
	        old[i][j] = smallmap[i-1][j-1];
	    }
    }

   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

   for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }

   /*
    *  Update for a fixed number of iterations
    */

  //maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;
  int anchange = 1;

  //while (step <= maxstep)
  while (anchange > 0)
    {
      nchange = 0;
      anchange = 0;
      MPI_Issend(&old[M][1],N,MPI_INT,right,tag1,comm2d,&request1);
      MPI_Recv(&old[0][1],N,MPI_INT,left,tag1,comm2d,&status1);
      MPI_Wait(&request1,&status1);

      MPI_Issend(&old[1][1],N,MPI_INT,left,tag2,comm2d,&request2);
      MPI_Recv(&old[M+1][1],N,MPI_INT,right,tag2,comm2d,&status2);
      MPI_Wait(&request2,&status2);

      MPI_Datatype hrowtype;
      MPI_Type_vector(M,1,N+2,MPI_INT,&hrowtype);
      MPI_Type_commit(&hrowtype);

      MPI_Issend(&old[1][N],1,hrowtype,up,tag3,comm2d,&request3);
      MPI_Recv(&old[1][0],1,hrowtype,down,tag3,comm2d,&status3);
      MPI_Wait(&request3,&status3);

      MPI_Issend(&old[1][1],1,hrowtype,down,tag4,comm2d,&request4);
      MPI_Recv(&old[1][N+1],1,hrowtype,up,tag4,comm2d,&status4);
      MPI_Wait(&request4,&status4);

      for (i=1; i<=M; i++)
	    {
	     for (j=1; j<=N; j++)
	     {
	      oldval = old[i][j];
	      newval = oldval;

	      /*
	       * Set new[i][j] to be the maximum value of old[i][j]
	       * and its four nearest neighbours
	       */

	      if (oldval != 0)
		    {
		      if (old[i-1][j] > newval) newval = old[i-1][j];
		      if (old[i+1][j] > newval) newval = old[i+1][j];
		      if (old[i][j-1] > newval) newval = old[i][j-1];
		      if (old[i][j+1] > newval) newval = old[i][j+1];

		      if (newval != oldval)
		      {
		        ++nchange;
		      }
		    }

	      new[i][j] = newval;
	     }
	    }

      /*
       *  Report progress every now and then
       */

      if (step % printfreq == 0)
	     {
	       printf("percolate: number of changes on step %d is %d\n",
		      step, nchange);
	     }
      MPI_Allreduce(&nchange,&anchange,1,MPI_INT,MPI_SUM,comm2d);
     // printf("anchange = %d\n ",anchange);
      /*
       *  Copy back in preparation for next step, omitting halos
       */

      for (i=1; i<=M; i++)
	     {
	       for (j=1; j<=N; j++)
	       {
	         old[i][j] = new[i][j];
	       }
	     }

      int smallmap_sum = 0;
      int map_sum;
      float map_average;
      for (i=1; i<=M; i++)
        {
          for (j=1; j<=N; j++)
          {
             smallmap_sum += old[i][j];
          }
        }
      MPI_Allreduce(&smallmap_sum,&map_sum,1,MPI_INT,MPI_SUM,comm2d);
      map_average = map_sum / (L * L);

      if( rank == 0) {
          if (step % printfreq == 0) {
              printf("percolate: the average of the map array on step %d "
                     "is %f\n", step, map_average);
          }
      }

      step++;
    } 

  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. maxstep
   *  is too small)
   */

  if (nchange != 0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
	     maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
	{
	  smallmap[i-1][j-1] = old[i][j];
	}
    }
  

  int tag = 0;
  MPI_Status status;
  //MPI_Gather(&smallmap[0][0], M*N, MPI_INT, &map[0][0], M*N, MPI_INT, 0, comm);
  //MPI_Ssend(&smallmap[0][0], M*N, MPI_INT, 0, tag, comm2d);
  /*
  MPI_Datatype rectype;
  MPI_Type_vector(M,N,N,MPI_INT,&rectype);
  MPI_Type_commit(&rectype);*/

  if(rank != 0){
     MPI_Ssend(&smallmap[0][0], M*N, MPI_INT, 0, tag, comm2d);
  }else{
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
