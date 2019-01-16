
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

#ifndef EXCHANGEHALO_HPP
#define EXCHANGEHALO_HPP

#ifndef HPCG_NO_MPI
#include <mpi.h>
#include "Geometry.hpp"
#include <cstdlib>
#include "SparseMatrix.hpp"
#include <type_traits>
#include "Vector.hpp"

/*!
  Communicates data that is at the border of the part of the domain assigned to this processor.

  @param[in]    A The known system matrix
  @param[inout] x On entry: the local vector entries followed by entries to be communicated; on exit: the vector with non-local entries updated by other processors
 */
template<class T>
void ExchangeHalo(const SparseMatrix & A, Vector<T> & x) {

  // Extract Matrix pieces

  local_int_t localNumberOfRows = A.localNumberOfRows;
  int num_neighbors = A.numberOfSendNeighbors;
  local_int_t * receiveLength = A.receiveLength;
  local_int_t * sendLength = A.sendLength;
  int * neighbors = A.neighbors;
  local_int_t totalToBeSent = A.totalToBeSent;
  local_int_t * elementsToSend = A.elementsToSend;

  int size, rank; // Number of MPI processes, My process ID
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //
  //  first post receives, these are immediate receives
  //  Do not wait for result to come, will do that at the
  //  wait call below.
  //

  int MPI_MY_TAG = 99;

  MPI_Request * request = new MPI_Request[num_neighbors];


  if (x.optimizationData && std::is_same<T, float>::value) {
    float * sendBuffer = (float*) A.sendBuffer;
    float * const xv = (float*)x.optimizationData;

    //
    // Externals are at end of locals
    //
    float * x_external = (float *) xv + localNumberOfRows;

    // Post receives first
    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_recv = receiveLength[i];
      MPI_Irecv(x_external, n_recv, MPI_FLOAT, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
      x_external += n_recv;
    }


    //
    // Fill up send buffer
    //

    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xv[elementsToSend[i]];

    //
    // Send to each neighbor
    //

    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_send = sendLength[i];
      MPI_Send(sendBuffer, n_send, MPI_FLOAT, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
      sendBuffer += n_send;
    }
  } else {
    double * sendBuffer = A.sendBuffer;
    double * const xv = x.values;

    //
    // Externals are at end of locals
    //
    double * x_external = (double *) xv + localNumberOfRows;

    // Post receives first
    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_recv = receiveLength[i];
      MPI_Irecv(x_external, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
      x_external += n_recv;
    }


    //
    // Fill up send buffer
    //

    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) sendBuffer[i] = xv[elementsToSend[i]];

    //
    // Send to each neighbor
    //

    // TODO: Thread this loop
    for (int i = 0; i < num_neighbors; i++) {
      local_int_t n_send = sendLength[i];
      MPI_Send(sendBuffer, n_send, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD);
      sendBuffer += n_send;
    }
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

  delete [] request;

  return;
}
#endif // HPCG_NO_MPI

#endif // EXCHANGEHALO_HPP
