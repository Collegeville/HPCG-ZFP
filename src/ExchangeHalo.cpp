
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
#include <cmath>
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

  double * const xv = x.values;

  //
  // Externals are at end of locals
  //
  double * x_external = (double *) xv + localNumberOfRows;


  //
  // Fill up send buffer
  //
  if (x.optimizationData) {
    // TODO: Thread this loop
    for (local_int_t i=0; i<totalToBeSent; i++) {
      double xBlock [BLOCK_SIZE];
      local_int_t src = elementsToSend[i];
      DecodeBlock(x, src/BLOCK_SIZE, xBlock);
      sendBuffer[i] = xBlock[src%BLOCK_SIZE];
    }
  } else {
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
  local_int_t total_recv = 0;
  for (int i = 0; i < num_neighbors; i++) {
    local_int_t n_recv = receiveLength[i];
    MPI_Irecv(x_external+total_recv, n_recv, MPI_DOUBLE, neighbors[i], MPI_MY_TAG, MPI_COMM_WORLD, request+i);
    total_recv += n_recv;
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
    double xBlock[BLOCK_SIZE];
    if (total_recv) {

      int offset = localNumberOfRows%BLOCK_SIZE;
      local_int_t block = localNumberOfRows/BLOCK_SIZE;
      DecodeBlock(x, block, xBlock);
      for (int j = offset; j < BLOCK_SIZE; j++) {
        xBlock[j] = x_external[j-offset];
      }
      EncodeBlock(x, block, xBlock);
      block++;
      local_int_t i = BLOCK_SIZE-offset;
      for (; block < ceil((total_recv+localNumberOfRows)/(double)BLOCK_SIZE); block++){
        EncodeBlock(x, block, x_external+i);
        i += BLOCK_SIZE;
      }
    }
  }

  delete [] request;

  return;
}
#endif
// ifndef HPCG_NO_MPI
