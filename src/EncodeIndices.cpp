
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

#include "EncodeIndices.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <set>
#include <queue>
#include <vector>
#include "CompressionData.hpp"


// computes floor(log_2(v))+1 where v is positive
inline int numBits(local_int_t v) {
  //copied from https://graphics.stanford.edu/~seander/bithacks.html
  const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
  const unsigned int S[] = {1, 2, 4, 8, 16};

  register unsigned int r = 0; // result of log2(v) will go here
  for (int i = 4; i >= 0; i--) {
    if (v & b[i]) {
      v >>= S[i];
      r |= S[i];
    }
  }
  //post loop r == floor(log2(v))
  return r+1;
}

typedef std::map<local_int_t, local_int_t> map_t;

struct HeapNode {
  bool terminal;
  local_int_t value;
  local_int_t count;

  HeapNode * leftChild;
  HeapNode * rightChild;

  HeapNode(local_int_t value, local_int_t count)
      : terminal(true), value(value), count(count),
        leftChild(0), rightChild(0) {}
  HeapNode(HeapNode * left, HeapNode * right)
      : terminal(false), value(0), leftChild(left), rightChild(right) {
    count = left->count + right->count;
  }

  ~HeapNode() {
    if (leftChild) {
      delete leftChild;
    }
    if (rightChild) {
       delete rightChild;
    }
  }
};


class compareHeapNodes {
   public:
   bool operator()(const HeapNode * node1, const HeapNode * node2)
   {
        return node1->count > node2->count;
   }
};

HeapNode* buildHeap(const map_t & indices) {
  std::priority_queue<HeapNode*, std::vector<HeapNode*>, compareHeapNodes> nodes;
  for (map_t::const_iterator it = indices.begin(); it != indices.end(); it++) {
    HeapNode * node = new HeapNode(it->first, it->second);
    nodes.push(node);
  }
  while (nodes.size() > 1) {
    HeapNode * node1 = nodes.top();
    nodes.pop();
    HeapNode * node2 = nodes.top();
    nodes.pop();
    nodes.push(new HeapNode(node1, node2));
  }
  return nodes.top();
}

void compressLevel(const HeapNode* root, std::vector<huffmanResult_t*> & tables,
                std::map<local_int_t, std::pair<uint8_t, uint64_t>> & encodings, const std::pair<uint8_t, uint64_t> prefix) {
  huffmanResult_t * currentTable = new huffmanResult_t[1<<WINDOW_SIZE];
  tables.push_back(currentTable);

  for (int i = (1 << WINDOW_SIZE)-1; i >= 0; i--) {
    const HeapNode * curNode = root;

    bool terminated = false;
    for (int j = 0; j < WINDOW_SIZE; j++) {
      if (curNode->terminal) {
        std::pair<uint8_t, uint64_t> encoding (prefix.first+j, prefix.second | ((i & ((1<<j)-1)) << prefix.first));
        encodings[curNode->value] = encoding;

        currentTable[i] = RESULT_DONE_MASK
                        | ((huffmanResult_t)j << RESULT_SIZE_OFFSET)
                        | curNode->value;
        terminated = true;
        break;
      } else {
        if (i & (1<<j)) {
          curNode = curNode->leftChild;
        } else {
          curNode = curNode->rightChild;
        }
      }
    }
    if (!terminated) {
      if (curNode->terminal) {
        std::pair<uint8_t, uint64_t> encoding (prefix.first+WINDOW_SIZE, prefix.second | ((i & ((1<<WINDOW_SIZE)-1)) << prefix.first));
        encodings[curNode->value] = encoding;

        currentTable[i] = RESULT_DONE_MASK
                        | ((huffmanResult_t)WINDOW_SIZE << RESULT_SIZE_OFFSET)
                        | curNode->value;
      } else {
        assert(tables.size() < std::numeric_limits<local_int_t>::max());
        currentTable[i] = tables.size();
        std::pair<uint8_t, uint64_t> newPrefix (prefix.first+WINDOW_SIZE, prefix.second | i<<prefix.first);
        compressLevel(curNode, tables, encodings, newPrefix);
      }
    }
  }
}


/*!rootPrefix
  Compresses a block into the given Vector.
  Note that the last block may be partial.

  @param[inout] vect   The vector to read from, must have optimization data.
  @param[in] blockid   The index of the block to read
  @param[in] block    The memory with the new content of the block.
*/
double EncodeIndices(SparseMatrix & mat) {

  double memUsage = 0;

  CompressionData & data = *(CompressionData*)mat.optimizationData;

  local_int_t maxFirst = 0;
  map_t offsetIndices;

  for (int i = 0; i < mat.localNumberOfRows; i++) {
    int nonzerosInRow = mat.nonzerosInRow[i];
    local_int_t * inds = mat.mtxIndL[i];
    double * vals = mat.matrixValues[i];
    double *& diagPtr = mat.matrixDiagonal[i];

    // need row indices in ascending order
    if (!std::is_sorted(inds, inds+nonzerosInRow)) {
      for (int j = 0; j<nonzerosInRow-1; j++) {
        local_int_t * nextElement =  //TODO build compressed tables & stream
        std::min_element(inds+j, inds+nonzerosInRow);
        int index = nextElement - inds;
        std::swap(vals[j], vals[index]);
        std::swap(inds[j], inds[index]);
        if (vals+index == diagPtr) {
          diagPtr = vals+j;
        } else if (vals+j == diagPtr) {
          diagPtr = vals+index;
        }
      }
    }

    local_int_t cur = inds[0];
    // NOTE: if element not present, it is initialized with it's default value (0)
    if (cur > maxFirst) {
      maxFirst = cur;
    }
    for (int j = 1; j < nonzerosInRow; j++) {
      local_int_t prev = cur;
      cur = inds[j];
      // NOTE: if element not present, it is initialized with it's default value (0)
      offsetIndices[cur-prev]++;
    }
  }

  data.firstIndBits = numBits(maxFirst);
  data.firstIndMask = (1 << data.firstIndBits)-1;

  HeapNode* offsetNodes = buildHeap(offsetIndices);
  std::vector<huffmanResult_t*> tableVector(0);
  std::map<local_int_t, std::pair<uint8_t, uint64_t>> encodings;
  compressLevel(offsetNodes, tableVector, encodings, std::pair<uint8_t, uint64_t>(0, 0));
  delete offsetNodes;

  data.numTables = tableVector.size();
  data.tables = new huffmanResult_t* [data.numTables];
  memUsage += sizeof(huffmanResult_t)*data.numTables;
  std::copy(tableVector.begin(), tableVector.end(), data.tables);

  uint8_t ** rows = new uint8_t* [mat.localNumberOfRows];
  memUsage += sizeof(uint8_t*)*mat.localNumberOfRows;
  for (local_int_t row = 0; row < mat.localNumberOfRows; row++) {
    const int tempRowSize = mat.nonzerosInRow[row]*8;
    uint8_t * tempRow = new uint8_t[tempRowSize](); //init to 0
    std::memcpy(tempRow, mat.mtxIndL[row], sizeof(local_int_t));
    int position = data.firstIndBits;

    local_int_t prevInd = mat.mtxIndL[row][0];
    for (local_int_t j = 1; j < mat.nonzerosInRow[row]; j++) {
      local_int_t curInd = mat.mtxIndL[row][j];
      local_int_t offset = curInd - prevInd;
      std::pair<uint8_t, uint64_t> encoding = encodings.at(offset);
      uint8_t bitCount = encoding.first;
      uint64_t bitValue = encoding.second;

      bitValue = (bitValue << position%8) | tempRow[position/8];
      std::memcpy(tempRow+position/8, &bitValue, sizeof(uint64_t));
      position += bitCount;

      prevInd = curInd;
    }

    //need to make sure the working array is large enough
    assert(position/8+8 < tempRowSize);

    rows[row] = new uint8_t[position/8+8];
    memUsage += sizeof(uint8_t)*(position/8+8);
    memcpy(rows[row], tempRow, position/8+8);
    delete [] tempRow;
  }
  data.compressed = rows;
}
