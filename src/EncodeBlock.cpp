
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
#include "Vector.hpp"


int withinTolerance(const double value, const double error) {
  if (errorMode & AND_MODES) {
    int withinTol = 1;
    if (errorMode & PW_REL) {
      withinTol &= error <= fabs(value)*PWREL_BOUND_RATIO;
    }
    if(errorMode & ABS) {
      withinTol &= error <= ABS_ERROR_BOUND;
    }
    return withinTol;
  } else {
    int withinTol = 0;
    if (errorMode & PW_REL) {
      withinTol |= error <= fabs(value)*PWREL_BOUND_RATIO;
    }
    if(errorMode & ABS) {
      withinTol |= error <= ABS_ERROR_BOUND;
    }
    return withinTol;
  }
}

inline void useUncompressed(const double value, double (& previousElements)[4], double * uncompressedArray, const int index, int & uncompressedCount){
  uncompressedArray[uncompressedCount] = value;
  previousElements[index%4] = value;
  uncompressedCount++;
}

inline void tryNeighborCompression(const double value, const double neighborErr, double (& previousElements)[4], char * compressed, double * uncompressedArray, const int index, int & uncompressedCount) {
  if (withinTolerance(value, neighborErr)) {
    compressed[index/4] |= NEIGHBOR << (2*(index%4));
    previousElements[index%4] = previousElements[(index-1)%4];
  } else {
    useUncompressed(value, previousElements, uncompressedArray, index, uncompressedCount);
  }
}

inline void tryLinearCompression(const double value, const double linPredict, const double linErr, double (& previousElements)[4], char * compressed, double * uncompressedArray,  const int index, int & uncompressedCount) {
  if (withinTolerance(value, linErr)) {
    compressed[index/4] |= LINEAR << (2*(index%4));
    previousElements[index%4] = linPredict;
  } else {
    useUncompressed(value, previousElements, uncompressedArray, index, uncompressedCount);
  }
}

inline void tryQuadraticCompression(const double value, const double quadPredict, const double quadErr, double (& previousElements)[4], char * compressed, double * uncompressedArray,  const int index, int & uncompressedCount) {
  if (withinTolerance(value, quadErr)) {
    compressed[index/4] |= QUADRATIC << (2*(index%4));
    previousElements[index%4] = quadPredict;
  } else {
    useUncompressed(value, previousElements, uncompressedArray, index, uncompressedCount);
  }
}

/*!
  Compresses a block into the given Vector.
  Note that the last block may be partial.

  @param[inout] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] block    The memory with the new content of the block.

  @return returns 0 upon success and 1 otherwise
*/
int EncodeBlock(Vector & vect, const local_int_t blockid, const double * block) {
  char * compressed = ((char *)vect.optimizationData)+BLOCK_BYTES*blockid;
  double * uncompressedArray = (double *)(compressed + COMPRESSED_BYTES);

  //length of the current block in values
  local_int_t blockLength = (blockid*BLOCK_SIZE <= vect.localLength) ? BLOCK_SIZE : vect.localLength%BLOCK_SIZE;

  double value, neighbor, neighborErr, linear, linErr, quad, quadErr;

  //start compressed section as all 0's
  for (int i = 0; i < COMPRESSED_BYTES; i++) {
    compressed[i] = 0;
  }

  //first value must be stored uncompressed
  uncompressedArray[0] = block[0];

  //Number of uncompressed values
  int uncompressedCount = 1;

  //need the last few elements to compress next element
  int index = 1;
  double previousElements [4] = {block[0], 0, 0, 0};

  value = block[index];
  //linear and quadratic doesn't work with 1 pnt
  tryNeighborCompression(value, fabs(block[0]-value), previousElements, compressed, uncompressedArray, index, uncompressedCount);
  index++;

  //quadratic doesn't work with 2 pnt
  value = block[index];
  linear = 2*previousElements[index-1] - previousElements[index-2];
  linErr = fabs(linear - value);
  neighbor = previousElements[index-1];
  neighborErr = fabs(neighbor - value);
  if (neighborErr <= linErr) {
    tryNeighborCompression(value, neighborErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
  } else {
    tryLinearCompression(value, linear, linErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
  }
  index++;


  for (; index < blockLength; index++) {
    value = block[index];
    linear = 2*previousElements[(index-1)%4] - previousElements[(index-2)%4];
    neighbor = previousElements[(index-1)%4];
    quad = 3*previousElements[(index-1)%4] - 3*previousElements[(index-2)%4] + previousElements[(index-3)%4];
    linErr = fabs(linear - value);
    neighborErr = fabs(neighbor - value);
    quadErr = fabs(quad - value);

    if (neighborErr <= linErr) {
      if(neighborErr <= quadErr){
        tryNeighborCompression(value, neighborErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
      } else {
        tryQuadraticCompression(value, quad, quadErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
      }
    } else {
      if(linErr <= quadErr){
        tryLinearCompression(value, linear, linErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
      } else {
        tryQuadraticCompression(value, quad, quadErr, previousElements, compressed, uncompressedArray, index, uncompressedCount);
      }
    }
  }

  return 0;
}
