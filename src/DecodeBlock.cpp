
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

#include "DecodeBlock.hpp"

#include <cassert>
#include "CompressionData.hpp"
#include "Vector.hpp"


//helper method with decode logic
inline int DoDecode(const local_int_t startIndex, const local_int_t blockLength, int & uncompressedCount, const char * compressed, double * block) {
  const double * uncompressedArray = (double *)(compressed + COMPRESSED_BYTES);

  for (int index = startIndex; index < blockLength; index++) {
    double value;
    switch((compressed[index/4] >> (2*(index&0b11)))&0b11) {
      case UNCOMPRESSED:
        value = uncompressedArray[uncompressedCount];
        uncompressedCount++;
        break;
      case NEIGHBOR:
        value = block[index-1];
        break;
      case LINEAR:
        value = 2*block[index-1] - block[index-2];
        break;
      case QUADRATIC:
        value = 3*block[index-1] - 3*block[index-2] + block[index-3];
        break;
      default:
        //This should never occur
        assert(1);
        return 1;
    }

    block[index] = value;
  }
  return 0;
}


/*!
  Decompresses a block from the given Vector.

  @param[in] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.

  @return returns 0 upon success and 1 otherwise
*/
int DecodeBlock(const Vector & vect, const local_int_t blockid, double * block) {
  const char * compressed = ((char *)vect.optimizationData)+BLOCK_BYTES*blockid;
  const double * uncompressedArray = (double *)(compressed + COMPRESSED_BYTES);

  //length of the current block in values
  local_int_t blockLength = (blockid*BLOCK_SIZE <= vect.localLength) ? BLOCK_SIZE : vect.localLength%BLOCK_SIZE;

  block[0] = uncompressedArray[0];
  int uncompressedCount = 1;

  DoDecode(1, blockLength, uncompressedCount, compressed, block);

  return 0;
}

/*!
  Decompresses part of a block from the given Vector.
  Only the first blockLength values are read.

  @param[in] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] blockLength
  @param[out] block    The memory to fill with the content of the block.

  @return returns 0 upon success and 1 otherwise
*/
int PartialDecodeBlock(const Vector & vect, const local_int_t blockid, const local_int_t blockLength, double * block) {
  const char * compressed = ((char *)vect.optimizationData)+BLOCK_BYTES*blockid;
  const double * uncompressedArray = (double *)(compressed + COMPRESSED_BYTES);

  block[0] = uncompressedArray[0];
  int uncompressedCount = 1;

  DoDecode(1, blockLength, uncompressedCount, compressed, block);

  if (blockLength != BLOCK_SIZE) {
    //If there is room at the end of the block, place the uncompressed count to allow resuming
    ((int*)(block+blockLength))[0] = uncompressedCount;
  }

  return 0;
}

/*!
  Resumes where PartialDecodeBlock left off.

  @param[in] vect         The vector to read from, must have optimization data.
  @param[in] blockid      The index of the block to read
  @param[in] blockLength
  @param[in] previousEnd  The length of the last decode
  @param[inout] block     The partially decoded block

  @return returns 0 upon success and 1 otherwise

*/
int ResumePartialDecodeBlock(const Vector & vect, const local_int_t blockid, const local_int_t blockLength, const local_int_t previousEnd, double * block) {
  const char * compressed = ((char *)vect.optimizationData)+BLOCK_BYTES*blockid;

  int uncompressedCount = ((int*)(block+previousEnd))[0];

  int ret = DoDecode(previousEnd, blockLength, uncompressedCount, compressed, block);

  if (blockLength != BLOCK_SIZE) {
    //If there is room at the end of the block, place the uncompressed count to allow resuming
    ((int*)(block+blockLength))[0] = uncompressedCount;
  }
  return ret;
}
