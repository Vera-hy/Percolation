#include <stdio.h>
#include <mpi.h>

#include "percolate.h"

/*
 *  Initialise MPI, compute the size
 */
void mp_start(int *size){

    MPI_Init(NULL, NULL);

    MPI_Comm_size(MPI_COMM_WORLD, size);
}

/*
 *  Finalize MPI
 */
int mp_stop(){
    MPI_Finalize();
    return 0;
}

/*
 * Blocking synchronous send.
 */
void mpSsend(void *sendbuf, int count, MPI_Datatype datatype, int dest,
        int tag, int comm){

    MPI_Ssend(sendbuf, count, datatype, dest, tag, comm);

}

/*
 * Determines process coords in cartesian topology given rank in group.
 */
void mpCartcoords(int comm, int rank, int *coord){

    MPI_Cart_coords(comm, rank, 2, coord);

}

/*
 * Blocking receive for a message.
 */
void mpRecv(void *recvbuf, int count, MPI_Datatype datatype, int source,
        int tag, int comm, MPI_Status *status){

    MPI_Recv(recvbuf, count, datatype, source, tag, comm, status);

}

/*
 * Broadcasts a message from the process with
 * rank "0" to all other processes of the communicator.
 */
void mpBcast(void *buffer, int count, int comm){

    MPI_Bcast(buffer, count, MPI_INT, 0, comm);

}

/*
 * Creates a division of processors in a cartesian grid.
 */
void mpDimscreate(int nnodes, int *dims){

    MPI_Dims_create(nnodes, ndims, dims);

}

/*
 * Makes a new communicator to which topology information has been attached.
 */
void mpCartcreate(int *dims, int *period, int reorder, int *comm){

    MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,period,reorder,comm);

}

/*
 * Determines the rank of the calling process in the communicator.
 */
void mpCommrank(int comm, int *rank){

    MPI_Comm_rank(comm,rank);

}

/*
 * Returns the shifted source and destination ranks,
 * given a shift direction and amount.
 */
void mpCartshift(int comm, int direction, int disp, int *source,
        int *dest ){

    MPI_Cart_shift(comm,direction,disp,source,dest);

}

/*
 * Starts a non-blocking synchronous send.
 */
void mpIssend(void *buf, int count, MPI_Datatype datatype, int dest, int tag,
        int comm, MPI_Request *request){

    MPI_Issend(buf, count, datatype, dest, tag, comm, request);

}

/*
 * Waits for all given MPI Requests to complete.
 */
void mpWaitall(  int count, MPI_Request array_of_requests[],
        MPI_Status array_of_statuses[]){

    MPI_Waitall(count, array_of_requests, array_of_statuses);

}

/*
 * Creates a vector datatype.
 */
void mpVector(int count, int blocklength,int stride, MPI_Datatype oldtype,
        MPI_Datatype *newtype){

    MPI_Type_vector(count, blocklength, stride, oldtype, newtype);

}

/*
 * Commits the datatype.
 */
void mpTypecommit(MPI_Datatype * datatype){

    MPI_Type_commit(datatype);

}
/*
 * Combines values from all processes and distributes the result
 * back to all processes.
 */
void mpgsum(void *sendbuf, void *recvbuf, int count, int comm){

    MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, comm);

}

/*
 * Starts a non-blocking receive.
 */
void mpIrecv(void *recvbuf, int count, MPI_Datatype datatype, int source,
        int tag, int comm, MPI_Request *request){

    MPI_Irecv(recvbuf, count, datatype, source, tag, comm, request);

}


