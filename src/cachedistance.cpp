// in cachedistance.cpp

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "sequence.h"
#include "cachedistance.h"

#include <limits.h>
#include <algorithm>
#include <memory>

// Boost symmetric matrix
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/foreach.hpp>

using namespace Rcpp;

CacheKey::CacheKey(const Sequence &s1, const Sequence &s2) : s1_(s1), s2_(s2) {
}

bool CacheKey::operator<(const CacheKey& k) const {
  if (s1_ < k.s1_)
    return true;
  else if (k.s1_ < s1_)
    return false;
  else
    return (s2_ < k.s2_);
}

CacheDistance::CacheDistance(const boost::shared_ptr< Distance > &dist) : dist_(dist) {

}

CacheDistance::~CacheDistance() {};

unsigned int CacheDistance::distance(const Sequence &sequence1, const Sequence &sequence2) {
  // FIXME: direction of distance can be optimized

  // Prevent the trivial case
  if (sequence1 == sequence2) {
    return 0;
  }

  bool have_return_val = false;
  unsigned int d = 0;

#pragma omp critical
  {
    // Find in Cache
    std::map< CacheKey, unsigned int >::iterator idx = cache_.find(CacheKey(sequence1, sequence2));
    if (idx != cache_.end())  {
      have_return_val = true;
      d = idx->second;
    }
  }

  if (have_return_val)
    return d;
  else {
    // Cache miss
    d = dist_->distance(sequence1, sequence2);
#pragma omp critical
    cache_.insert(std::pair< CacheKey, unsigned int>(CacheKey(sequence1, sequence2), d));
  }

  return(d);
}

unsigned int CacheDistance::min_seq_distance(const std::vector<Sequence> &seqs, const Sequence &sequence, const size_t n) {
  size_t n_elements = seqs.size();

  unsigned int global_min_dist = UINT_MAX;

  for (size_t seq1_index = 0; seq1_index < n_elements; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance  = distance(sequence1,sequence);

    if (min_distance < global_min_dist)
      global_min_dist = min_distance;

  }

  return(global_min_dist);
}

unsigned int CacheDistance::min_set_distance(const std::vector<Sequence> &sequences, const size_t n) {
  size_t n_elements = sequences.size();

  unsigned int set_min_distance = UINT_MAX;

  for (size_t row_idx(0); row_idx < n_elements; row_idx++) {
    for (size_t col_idx(row_idx+1); col_idx < n_elements; col_idx++) {

      Sequence sequence1 = sequences.at(row_idx);
      Sequence sequence2 = sequences.at(col_idx);

      unsigned int min_distance  = distance(sequence1, sequence2);

      if (min_distance < set_min_distance)
        set_min_distance = min_distance;
    }
  }
  return(set_min_distance);
}

bool CacheDistance::is_seq_insertable(const std::vector<Sequence> &seqs, const Sequence &sequence, const size_t n, const unsigned int min_dist) {
  size_t n_elements = seqs.size();

  bool seq_is_insertable = true;

  for (size_t seq1_index = 0; seq1_index < n_elements && seq_is_insertable; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance  = distance(sequence1, sequence);

    if (min_distance < min_dist)
      seq_is_insertable = false;
  }

  return seq_is_insertable;
}

Rcpp::DataFrame CacheDistance::demultiplex(const std::vector< std::string > &barcodes, const std::vector< std::string > &reads) {
  return(dist_->demultiplex(barcodes, reads));
}

