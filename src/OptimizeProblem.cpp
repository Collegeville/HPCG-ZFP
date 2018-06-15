
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
 @file OptimizeProblem.cpp

 HPCG routine
 */

#include "OptimizeProblem.hpp"
#include <cassert>
#include <cmath>
#include <cstdint>
#include "CompressionData.hpp"
#include <iostream>
#include "zfparray1.h"

double optimization_allocation = 0.0;


/*!
  helper function to create ZFP arrays for vectors

  @param[inout] vect   The vector to create a zfp array for

  @return returns the number of bytes used for the array
*/
int CreateZFPArray(SparseMatrix & mat){

  CompressionData * data = new CompressionData();
  mat.optimizationData = data;

  data->bufferSize = blockRate * ceil(mat.localNumberOfNonzeros/4.0);
  data->buffer = new uint8_t[data->bufferSize];
  assert(data->buffer);
  data->stream = stream_open(data->buffer, data->bufferSize);
  assert(data->stream);

  zfp_stream * zfp = zfp_stream_open(data->stream);
  assert(zfp);
  data->zfp = zfp;
  double actualRate = zfp_stream_set_rate(zfp, bitRate, zfp_type_double, 1, 1);
  assert(actualRate == (double)bitRate);

  local_int_t * rowStarts = new local_int_t[mat.localNumberOfRows];
  data->rowStarts = rowStarts;
  local_int_t * diagonalIndices = new local_int_t[mat.localNumberOfRows];
  data->diagonalIndices = diagonalIndices;

  local_int_t position = 0;
  double currentBlock [4];
  for (int i = 0; i < mat.localNumberOfRows; i++) {
    int nnz = mat.nonzerosInRow[i];
    rowStarts[i] = position;
    diagonalIndices[i] = rowStarts[i] + (mat.matrixDiagonal[i] - mat.matrixValues[i]);
    for (int j = 0; j < nnz; j++) {
      currentBlock[(position+j)%4] = mat.matrixValues[i][j];
      if ((position+j)%4 == 3) {
        EncodeFullBlock(mat, (position+j)/4, currentBlock);
      }
    }
    position += nnz;
  }
  if (position%4) {
    EncodeBlock(mat, position/4, currentBlock);
  }

  return sizeof(*data)
          + data->bufferSize //buffer
          + sizeof(uint)+stream_word_bits+3*sizeof(uint64*); //sizeof(*data->stream)
          + sizeof(*data->zfp)
          + sizeof(local_int_t)*2*mat.localNumberOfRows; //rowStarts & diagonalIndices
}

/*!
  Optimizes the data structures used for CG iteration to increase the
  performance of the benchmark version of the preconditioned CG algorithm.

  @param[inout] A      The known system matrix, also contains the MG hierarchy in attributes Ac and mgData.
  @param[inout] data   The data structure with all necessary CG vectors preallocated
  @param[inout] b      The known right hand side vector
  @param[inout] x      The solution vector to be computed in future CG iteration
  @param[inout] xexact The exact solution vector

  @return returns 0 upon success and non-zero otherwise

  @see GenerateGeometry
  @see GenerateProblem
*/
int OptimizeProblem(SparseMatrix & A, CGData & data, Vector & b, Vector & x, Vector & xexact) {

  // This function can be used to completely transform any part of the data structures.
  // Right now it does nothing, so compiling with a check for unused variables results in complaints

#if defined(HPCG_USE_MULTICOLORING)
  const local_int_t nrow = A.localNumberOfRows;
  std::vector<local_int_t> colors(nrow, nrow); // value `nrow' means `uninitialized'; initialized colors go from 0 to nrow-1
  int totalColors = 1;
  colors[0] = 0; // first point gets color 0

  // Finds colors in a greedy (a likely non-optimal) fashion.

  for (local_int_t i=1; i < nrow; ++i) {
    if (colors[i] == nrow) { // if color not assigned
      std::vector<int> assigned(totalColors, 0);
      int currentlyAssigned = 0;
      const local_int_t * const currentColIndices = A.mtxIndL[i];
      const int currentNumberOfNonzeros = A.nonzerosInRow[i];

      for (int j=0; j< currentNumberOfNonzeros; j++) { // scan neighbors
        local_int_t curCol = currentColIndices[j];
        if (curCol < i) { // if this point has an assigned color (points beyond `i' are unassigned)
          if (assigned[colors[curCol]] == 0)
            currentlyAssigned += 1;
          assigned[colors[curCol]] = 1; // this color has been used before by `curCol' point
        } // else // could take advantage of indices being sorted
      }

      if (currentlyAssigned < totalColors) { // if there is at least one color left to use
        for (int j=0; j < totalColors; ++j)  // try all current colors
          if (assigned[j] == 0) { // if no neighbor with this color
            colors[i] = j;
            break;
          }
      } else {
        if (colors[i] == nrow) {
          colors[i] = totalColors;
          totalColors += 1;
        }
      }
    }
  }

  std::vector<local_int_t> counters(totalColors);
  for (local_int_t i=0; i<nrow; ++i)
    counters[colors[i]]++;

  local_int_t old, old0;
  for (int i=1; i < totalColors; ++i) {
    old0 = counters[i];
    counters[i] = counters[i-1] + old;
    old = old0;
  }
  counters[0] = 0;

  // translate `colors' into a permutation
  for (local_int_t i=0; i<nrow; ++i) // for each color `c'
    colors[i] = counters[colors[i]]++;
#endif

  double bytes_used = 0;

  SparseMatrix * Anext = &A;

  while (Anext) {
    bytes_used += CreateZFPArray(*Anext);
    Anext = Anext->Ac;
  }
  optimization_allocation = bytes_used;

  return 0;
}

// Helper function (see OptimizeProblem.hpp for details)
double OptimizeProblemMemoryUse(const SparseMatrix & A) {

  return optimization_allocation;

}
