
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
 @file ComputeProlongation.cpp

 HPCG routine
 */

#include "ComputeProlongation.hpp"

/*!
  Routine to compute the coarse residual vector.

  @param[in]  Af - Fine grid sparse matrix object containing pointers to current coarse grid correction and the f2c operator.
  @param[inout] xf - Fine grid solution vector, update with coarse grid correction.

  Note that the fine grid residual is never explicitly constructed.
  We only compute it for the fine grid points that will be injected into corresponding coarse grid points.

  @return Returns zero on success and a non-zero value otherwise.
*/
int ComputeProlongation(const SparseMatrix & Af, Vector & xf) {

  Vector & xc = *Af.mgData->xc;
  const local_int_t * f2c = Af.mgData->f2cOperator;
  const local_int_t nc = Af.mgData->rc->localLength;

  double cBlock[BLOCK_SIZE];
  double fBlock[BLOCK_SIZE];
  local_int_t currentFineBlock = f2c[0]/BLOCK_SIZE;
  DecodeBlock(xf, currentFineBlock, fBlock);

  for (local_int_t block = 0; block < nc/BLOCK_SIZE; block++){
    local_int_t i = block*BLOCK_SIZE;
    PartialDecodeBlock(xc, block, BLOCK_SIZE, cBlock);
    for (local_int_t j = 0; j < BLOCK_SIZE; j++) {
      local_int_t dest = f2c[i+j];
      if (dest/BLOCK_SIZE != currentFineBlock) {
        //need to encode previous block
        EncodeBlock(xf, currentFineBlock, fBlock);
        currentFineBlock = dest/BLOCK_SIZE;
        DecodeBlock(xf, currentFineBlock, fBlock);
      }
      fBlock[dest%BLOCK_SIZE] = cBlock[j];
    }
  }
  if (nc%BLOCK_SIZE) {
    local_int_t block = nc/BLOCK_SIZE;
    local_int_t i = block*BLOCK_SIZE;
    DecodeBlock(xc, block, cBlock);
    for (local_int_t j = 0; j < nc%BLOCK_SIZE; j++) {
      local_int_t dest = f2c[i+j];
      if (dest/BLOCK_SIZE != currentFineBlock) {
        //need to encode previous block before decoding next one
        EncodeBlock(xf, currentFineBlock, fBlock);
        currentFineBlock = dest/BLOCK_SIZE;
        DecodeBlock(xf, currentFineBlock, fBlock);
      }
      fBlock[dest%BLOCK_SIZE] = cBlock[j];
    }
  }

  //Encode last block needed
  EncodeBlock(xf, currentFineBlock, fBlock);


  return 0;
}
