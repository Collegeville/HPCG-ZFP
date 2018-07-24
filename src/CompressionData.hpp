
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
 @file CGData.hpp

 HPCG data structure
 */

#ifndef COMPRESSIONDATA_HPP
#define COMPRESSIONDATA_HPP

#include <cmath>
#include <cstdint>
#include <cstring>
#include "Geometry.hpp"


constexpr int WINDOW_SIZE = 8;
constexpr int WINDOW_SIZE_BITS = 4; //bytes needed to hold a window size value
static_assert(std::pow(2, WINDOW_SIZE_BITS)-1 >= WINDOW_SIZE, "Not enought bits for window size");


typedef uint64_t huffmanResult_t;

constexpr huffmanResult_t HUFFMAN_RESULT_BITS = 8*sizeof(huffmanResult_t);
static_assert(8*sizeof(local_int_t) + WINDOW_SIZE_BITS + 1 <= HUFFMAN_RESULT_BITS, "Need room for flags and stuff");

constexpr huffmanResult_t RESULT_DONE_MASK = ((huffmanResult_t)1 << (HUFFMAN_RESULT_BITS - 1));
constexpr huffmanResult_t RESULT_SIZE_OFFSET = (HUFFMAN_RESULT_BITS-1-WINDOW_SIZE_BITS);
constexpr huffmanResult_t RESULT_SIZE_MASK = (((huffmanResult_t)1 << WINDOW_SIZE_BITS)-1) << RESULT_SIZE_OFFSET;
constexpr huffmanResult_t RESULT_VALUE_MASK = ((huffmanResult_t)1 << (8*sizeof(local_int_t)))-1;


constexpr int WINDOW_MASK = (1 << WINDOW_SIZE) -1;


struct CompressionData {
  huffmanResult_t ** tables;
  local_int_t numTables;
  uint8_t ** compressed;
  local_int_t numRows;

  uint64_t firstIndMask;
  int firstIndBits;
};


inline void DeleteCompressionData(CompressionData & data){
  if (data.tables) {
    for (local_int_t i = 0; i < data.numTables; i++) {
      delete [] data.tables[i];
    }
    delete [] data.tables;
    data.tables = 0;
  }

  if (data.compressed) {
    for (local_int_t i = 0; i < data.numRows; i++) {
      delete [] data.compressed[i];
    }
    delete [] data.compressed;
  }
}


inline huffmanResult_t peek(const huffmanResult_t * huffmanTable, const uint8_t * compressed, const uint64_t position) {
  uint32_t window;
  std::memcpy(&window, compressed+position/8, sizeof(window));
  window >>= position%8;

  return huffmanTable[window & WINDOW_MASK];
}

#endif // COMPRESSIONDATA_HPP
