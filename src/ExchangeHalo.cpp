
//@HEADER
// ***************************************************
//
// HPCG: High Performance Conjugate Gradient Benchmark
//
// Contact:
// Michael A. Heroux ( maherou@sandia.gov)
// Jack Dongarra     (dongarra@eecs.utk.edu)
// Piotr Luszczek    (luszczek@eecs.utk.edu)
//
// ***************************************************
//@HEADER

/*!
 @file ExchangeHalo.cpp

 HPCG routine
 */

// Compile this routine only if running with MPI
#ifndef HPCG_NO_MPI
#include <mpi.h>
#include "Geometry.hpp"
#include "ExchangeHalo.hpp"
#include <cstdlib>

/*!
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated; on exit: the vector with non-local entries updated by other processors
 */
void ExchangeHalo(const SparseMatrix & A, Vector & x) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  int num_neighbors = A.numberOfSendNeighbors;
  local_int_t * receiveLength = A.receiveLength;
  local_int_t * sendLength = A.sendLength;
  int * neighbors = A.neighbors;
  double * sendBuffer = A.sendBuffer;
  local_int_t totalToBeSent = A.totalToBeSent;
  local_int_t * elementsToSend = A.elementsToSend;

  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  double * x_external;

  if (x.optimizationData){
    zfp::array1d & xopt = *(zfp::array1d*)x.optimizationData;

    x_external = new double[x.localLength - localNumberOfRows];

    //
    // Fill up send buffer
    //

    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xopt[elementsToSend[i]];
  } else {
    double * const xv = x.values;

    //
    // Externals are at end of locals
    //
    x_external = (double *) xv + localNumberOfRows;

    //
    // Fill up send buffer
    //

    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xv[elementsToSend[i]];
  }


  //
  //  first post receives, these are immediate receives
  //  Do not wait for result to come, will do that at the
  //  wait call below.
  //

  int MPI_MY_TAG = 99;

  MPI_Request * request = new MPI_Request[num_neighbors];

  // Post receives first
  // TODO: Thread this loop
  double * x_recv_location = x_external;
  for (int i = 0; i < num_neighbors; i++) {
    local_int_t n_recv = receiveLength[i];
    MPI_Irecv(x_recv_location, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
    x_recv_location += n_recv;
  }

  //
  // Send to each neighbor
  //

  // TODO: Thread this loop
  for (int i = 0; i < num_neighbors; i++) {
    local_int_t n_send = sendLength[i];
    MPI_Send(sendBuffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
    sendBuffer += n_send;
  }

  //
  // Complete the reads issued above
  //

  MPI_Status status;
  // TODO: Thread this loop
  for (int i = 0; i < num_neighbors; i++) {
    if ( MPI_Wait(request+i, &status) ) {
      std::exit(-1); // TODO: have better error exit
    }
  }

  if (x.optimizationData){
    zfp::array1d & xopt = *(zfp::array1d*)x.optimizationData;
    for (int i = 0; i<x_recv_location - x_external; i++){
      xopt[i+localNumberOfRows] = x_external[i];
    }
    delete [] x_external;
  }

  delete [] request;

  return;
}
#endif
// ifndef HPCG_NO_MPI
