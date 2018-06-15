
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

#include "CompressionData.hpp"
#include "SparseMatrix.hpp"


/*!
  Decompresses a block from the given Vector.
  Note that the last block may be partial.

  @param[in] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.  Must have space for 4 doubles.
  @param[out] count    The location is set to number of values in the block

  @return returns the count, or 0 on an error
*/
int DecodeBlock(const SparseMatrix & mat, const local_int_t blockid, double * block) {
  CompressionData * data = (CompressionData*)mat.optimizationData;

  return DecodeBlock(data->zfp, mat.localNumberOfNonzeros, blockid, block);
}


/*!
  Decompresses a block from the given zfp_stream.

  @param[in] zfp       The zfp_stream to read from
  @param[in] arrayLength  The length of the data, used to determine if the requested block is a full or partial block
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.  Must have space for 4 doubles.

  @return returns 1 upon success and 0 otherwise

  @see DecodeBlock
*/
int DecodeBlock(zfp_stream * zfp, local_int_t arrayLength, local_int_t blockid, double * block) {
  bitstream * stream = zfp_stream_bit_stream(zfp);
  stream_rseek(stream, blockid * blockRate*8);

  int count = 4;
  //last block might be partial
  if (blockid-1 == arrayLength/4 && arrayLength%4 != 0) {
    count = arrayLength%4;
    if (zfp_decode_partial_block_strided_double_1(zfp, block, count, 1)) {
      count = 0;
    }
  } else {
    if (zfp_decode_block_double_1(zfp, block)) {
      count = 0;
    }
  }
  return count;
}


/*!
  Decompresses a full block from the given Vector.

  @param[in] vect      The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.  Must have space for 4 doubles.

  @return returns 1 upon success and 0 otherwise

  @see DecodeBlock
*/
int DecodeFullBlock(const SparseMatrix & mat, const local_int_t blockid, double * block) {
  CompressionData * data = (CompressionData*)mat.optimizationData;

  return DecodeFullBlock(data->zfp, blockid, block);
}

/*!
  Decompresses a full block from the given zfp_stream.

  @param[in] zfp       The zfp_stream to read from
  @param[in] blockid   The index of the block to read
  @param[out] block    The memory to fill with the content of the block.  Must have space for 4 doubles.

  @return returns 1 upon success and 0 otherwise

  @see DecodeBlock
*/
int DecodeFullBlock(zfp_stream * zfp, const local_int_t blockid, double * block) {
  bitstream * stream = zfp_stream_bit_stream(zfp);
  stream_rseek(stream, blockid * blockRate*8);

  return !zfp_decode_block_double_1(zfp, block);
}
