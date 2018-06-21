
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
#include "CompressionData.hpp"
#include "EncodeValues.hpp"

double optimizationAllocation = 0.0;


/*!
  helper function to create compressionData

  @param[inout] mat   The matrix to optimize

  @return the number of bytes used
*/
int CreateCompressedArray(SparseMatrix & mat){

  local_int_t neededCompression = ceil(mat.localNumberOfNonzeros/(double)VALUES_PER_COMPRESSED_BYTE);

  CompressionData * data = new CompressionData();
  mat.optimizationData = data;

  data->fValsCompressed = new unsigned char[neededCompression];
  assert(data->fValsCompressed);
  data->fValsUncompressed = new double[mat.localNumberOfNonzeros];
  assert(data->fValsUncompressed);
  data->bValsCompressed = new unsigned char[neededCompression];
  assert(data->bValsCompressed);
  data->bValsUncompressed = new double[mat.localNumberOfNonzeros];
  assert(data->bValsUncompressed);

  //reuse matrix diagonal allocation
  assert(sizeof(double**) >= sizeof(double)); //ensure the matrixDiagonal array is larger than our new copy.
  data->diagonalValues = (double*)mat.matrixDiagonal; //new double[mat.localNumberOfRows];
  assert(data->diagonalValues);

  EncodeValues(mat, true);
  EncodeValues(mat, false);

  //copy diagonal into array to get 8 vals per cache line instead of 1 val per cache line of mat.matrixDiagonal
  //The values are compressed serially, so a pointer/index of an entry doesn't work.
  for (local_int_t i = 0; i < mat.localNumberOfRows; i++){
    data->diagonalValues[i] = mat.matrixDiagonal[i][0];
  }


  return sizeof(*data)  //structure itself
          + 2*sizeof(unsigned char)*neededCompression //compressed arrays
          + 2*sizeof(local_int_t)*mat.localNumberOfNonzeros // uncompressed arrays
          + sizeof(double)*mat.localNumberOfRows; //diagonal
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

  optimizationAllocation = 0;

  SparseMatrix * Anext = &A;
  while (Anext) {
    optimizationAllocation += CreateCompressedArray(*Anext);
    Anext = Anext->Ac;
  }

  return 0;
}

// Helper function (see OptimizeProblem.hpp for details)
double OptimizeProblemMemoryUse(const SparseMatrix & A) {

  return optimizationAllocation;

}
