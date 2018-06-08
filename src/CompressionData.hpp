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

#include <cmath>

#ifndef COMPRESSIONDATA_HPP
#define COMPRESSIONDATA_HPP
/*!
  The number of values stored per block
*/
constexpr const int BLOCK_SIZE = 16;

/*!
  Whether to cache align the starts of blocks
*/
constexpr const int CACHE_ALIGNED = 0;

/*!
  The size of a cache line in bytes
*/
constexpr const int CACHE_LINE_SIZE = 64;

/**
  The number of bytes needed for the compressed part of each block
*/
constexpr const int COMPRESSED_BYTES = ceil(BLOCK_SIZE/4.0);

//used to calculate blockBytes
constexpr const int MIN_BLOCK_BYTES = COMPRESSED_BYTES + 8*BLOCK_SIZE;

/*!
  The number of bytes per block.
*/
constexpr const int BLOCK_BYTES = CACHE_ALIGNED ? ceil(MIN_BLOCK_BYTES/CACHE_LINE_SIZE)*CACHE_LINE_SIZE : MIN_BLOCK_BYTES;



//error bound modes
constexpr const int PW_REL = 0b0001;
constexpr const int ABS = 0b0010;
//constexpr const int REL = 0b0100;
constexpr const int AND_MODES = 0b1000;
constexpr const int ABS_AND_PW_REL = AND_MODES | ABS | PW_REL;
constexpr const int ABS_OR_PW_REL = AND_MODES | ABS | PW_REL;


/*!
  The error bound mode as per SZ
*/
constexpr const int errorMode = PW_REL;

//error bounds as per SZ
constexpr const double ABS_ERROR_BOUND = 1e-6;
constexpr const double PWREL_BOUND_RATIO = 1e-6;


constexpr const int UNCOMPRESSED = 0x00;
constexpr const int NEIGHBOR = 0x01;
constexpr const int LINEAR = 0x02;
constexpr const int QUADRATIC = 0x03;


#endif // COMPRESSIONDATA_HPP
