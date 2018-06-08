
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


/*!
  Decompresses a block from the given Vector.

  @param[in] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.

  @return returns 0 upon success and 1 otherwise
*/
int DecodeBlock(const Vector & vect, const local_int_t blockid, double * block) {

  //length of the current block in values
  local_int_t blockLength = (blockid*BLOCK_SIZE <= vect.localLength) ? BLOCK_SIZE : vect.localLength%BLOCK_SIZE;

  return PartialDecodeBlock(vect, blockid, blockLength, block);
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

  //Number of uncompressed values
  int uncompressedCount = 1;

  double previousElements [4] = {block[0], 0, 0, 0};

  for (int index = 1; index < blockLength; index++) {
    double value;
    switch((compressed[index/4] >> (2*(3-index%4)))&0x03) {
      case UNCOMPRESSED:
        value = uncompressedArray[uncompressedCount];
        uncompressedCount++;
        break;
      case NEIGHBOR:
        value = previousElements[(index-1)%4];
        break;
      case LINEAR:
        value = 2*previousElements[(index-1)%4] - previousElements[(index-2)%4];
        break;
      case QUADRATIC:
        value = 3*previousElements[(index-1)%4] - 3*previousElements[(index-2)%4] + previousElements[(index-3)%4];
        break;
      default:
        //This should never occur
        assert(1);
        return 1;
    }

    previousElements[index%4] = value;
    block[index] = value;
  }

  return 0;
}
