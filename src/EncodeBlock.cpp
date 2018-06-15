
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

#include "EncodeBlock.hpp"

#include "CompressionData.hpp"
#include "SparseMatrix.hpp"


/*!
  Compresses a block into the given Vector.
  Note that the last block may be partial.

  @param[inout] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] block    The memory with the new content of the block.

  @return returns 1 upon success and 0 otherwise
*/
int EncodeBlock(SparseMatrix & mat, const local_int_t blockid, double * block) {
  CompressionData * data = (CompressionData*)mat.optimizationData;

  return EncodeBlock(data->zfp, mat.localNumberOfNonzeros, blockid, block);
}


/*!
  Compresses a block into the given zfp_stream.
  Note that the last block may be partial.

  @param[in] zfp         The zfp_stream to read from, must have optimization data.
  @param[in] arrayLength The length of the data
  @param[in] blockid     The index of the block to read
  @param[in] block       The memory with the new content of the block.

  @return returns 1 upon success and 0 otherwise
*/
int EncodeBlock(zfp_stream * zfp, local_int_t arrayLength, local_int_t blockid, double * block) {
  bitstream * stream = zfp_stream_bit_stream(zfp);
  stream_wseek(stream, blockid * blockRate*8);

  int err;
  //last block might be partial
  if (blockid-1 == arrayLength/4 && arrayLength%4 != 0) {
    int count = arrayLength%4;
    err = !zfp_encode_partial_block_strided_double_1(zfp, block, count, 1);
  } else {
    err = !zfp_encode_block_double_1(zfp, block);
  }
  return err;
}


/*!
  Compresses a full block into the given Vector.

  @param[inout] vect   The vector to write to, must have optimization data.
  @param[in] blockid   The index of the block to write
  @param[in] block    The memory with the new content of the block.  Must contain 4 doubles.

  @return returns 1 upon success and 0 otherwise

  @see EncodeBlock
*/
int EncodeFullBlock(SparseMatrix & mat, const local_int_t blockid, double * block) {
  CompressionData * data = (CompressionData*)mat.optimizationData;

  return EncodeFullBlock(data->zfp, blockid, block);
}


/*!
  Compresses a full block into the given zfp_stream.
  Note that the last block may be partial.

  @param[in] zfp         The zfp_stream to read from, must have optimization data.
  @param[in] blockid     The index of the block to read
  @param[in] block       The memory with the new content of the block.

  @return returns 1 upon success and 0 otherwise
*/
int EncodeFullBlock(zfp_stream * zfp, local_int_t blockid, double * block) {
  bitstream * stream = zfp_stream_bit_stream(zfp);
  stream_wseek(stream, blockid * blockRate*8);

  return !zfp_encode_block_double_1(zfp, block);
}
