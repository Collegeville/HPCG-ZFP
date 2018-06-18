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
#include "Geometry.hpp"

#ifndef COMPRESSIONDATA_HPP
#define COMPRESSIONDATA_HPP

struct CompressionData_STRUCT {
  unsigned char * forwardCompressed;
  local_int_t * forwardUncompressed;
  unsigned char * backwardCompressed;
  local_int_t * backwardUncompressed;
};
typedef CompressionData_STRUCT CompressionData;

inline void DeleteCompressionData(CompressionData & data) {
  if (data.forwardCompressed) {
    delete data.forwardCompressed;
    data.forwardCompressed = 0;
  }
  if (data.forwardUncompressed) {
    delete data.forwardUncompressed;
    data.forwardUncompressed = 0;
  }
  if (data.backwardCompressed) {
    delete data.backwardCompressed;
    data.backwardCompressed = 0;
  }
  if (data.backwardUncompressed) {
    delete data.backwardUncompressed;
    data.backwardUncompressed = 0;
  }
}

constexpr const int UNCOMPRESSED = 0x00;
constexpr const int INCREMENT = 0x01;


#endif // COMPRESSIONDATA_HPP
