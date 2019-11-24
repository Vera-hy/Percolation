
#ifndef MPLIB_H
#define MPLIB_H

/*
 *  Prototypes for functions of MPI library.
 */

/*
 *  Initialise MPI and compute the size
 */
void mp_start(int *size);

/*
 *  Finalize MPI
 */
int mp_stop();

/*
 * Blocking synchronous send.
 */
void mpSsend(void *sendbuf, int count, MPI_Datatype datatype, int dest,
             int tag, int comm);

/*
 * Determines process coords in cartesian topology given rank in group.
 */
void mpCartcoords(int comm, int rank, int *coord);

/*
 * Blocking receive for a message.
 */
void mpRecv(void *recvbuf, int count, MPI_Datatype datatype, int source,
            int tag, int comm, MPI_Status *status);

/*
 * Broadcasts a message from the process with
 * rank "0" to all other processes of the communicator.
 */
void mpBcast(void *buffer, int count, int comm);

/*
 * Creates a division of processors in a cartesian grid.
 */
void mpDimscreate(int nnodes, int *dims);

/*
 * Makes a new communicator to which topology information has been attached.
 */
void mpCartcreate(int *dims, int *period, int reorder, int *comm);

/*
 * Determines the rank of the calling process in the communicator.
 */
void mpCommrank(int comm, int *rank);

/*
 * Returns the shifted source and destination ranks,
 * given a shift direction and amount.
 */
void mpCartshift(int comm, int direction, int disp, int *source,
                 int *dest );

/*
 * Starts a non-blocking synchronous send.
 */
void mpIssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
              int comm, MPI_Request *request);

/*
 * Waits for all given MPI Requests to complete.
 */
void mpWaitall(  int count, MPI_Request array_of_requests[],
                 MPI_Status array_of_statuses[]);

/*
 * Creates a vector datatype
 */
void mpVector(int count, int blocklength,int stride, MPI_Datatype oldtype,
              MPI_Datatype *newtype);

/*
 * Commits the datatype
 */
void mpTypecommit(MPI_Datatype * datatype);

/*
 * Combines values from all processes and distributes the result
 * back to all processes.
 */
void mpgsum(void *sendbuf, void *recvbuf, int count, int comm);

/*
 * Starts a non-blocking receive.
 */
void mpIrecv(void *recvbuf, int count, MPI_Datatype datatype, int source,
             int tag, int comm, MPI_Request *request);

/*
 * Determines process rank in communicator given Cartesian location.
 */
void mpCartrank(int comm, int *coords, int *rank);

#endif
