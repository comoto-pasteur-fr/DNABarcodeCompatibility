// in hammingdistance.cpp

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "sequence.h"
#include "hammingdistance.h"

#include <limits.h>
#include <algorithm>
#include <memory>

using namespace Rcpp;

HammingDistance::HammingDistance(const unsigned int cost_sub) : cost_sub_(cost_sub)  {
}

HammingDistance::~HammingDistance() {};

unsigned int HammingDistance::distance(const Sequence &seq1, const Sequence &seq2) {
    return(static_distance(seq1, seq2, cost_sub_));
}

unsigned int HammingDistance::min_seq_distance(const std::vector< Sequence > &seqs, const Sequence &sequence, const size_t n) {
    return(static_min_seq_distance( seqs, sequence, n, cost_sub_));
}

unsigned int HammingDistance::min_set_distance(const std::vector< Sequence > &sequences, const size_t n) {
    return(static_min_set_distance(sequences, n, cost_sub_));
}

bool HammingDistance::is_seq_insertable(const std::vector< Sequence > &seqs, const Sequence &sequence, const size_t n, const unsigned int min_dist) {
  size_t n_elements = seqs.size();

  bool seq_is_insertable = true;

  for (size_t seq1_index = 0; seq1_index < n_elements && seq_is_insertable; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance  = static_distance(sequence1, sequence, cost_sub_);

    if (min_distance < min_dist)
      seq_is_insertable = false;
  }

  return seq_is_insertable;

}

Rcpp::DataFrame HammingDistance::demultiplex(const std::vector< std::string > &barcodes, const std::vector< std::string > &reads) {
  if ((barcodes.size() < 2) || (reads.size() < 1)) {
    Rcpp::stop("At least one read and two barcodes need to be provided");
  }

  size_t n = barcodes[1].length();

  for (size_t i = 1; i < barcodes.size(); i++) {
    if (barcodes[i].length() != n) {
      Rcpp::stop("Length of all barcodes must be equal.");
    }
  }
  for (size_t i = 0; i < reads.size(); i++) {
    if (reads[i].length() != n) {
      Rcpp::stop("Length of all reads must be equal to barcode length.");
    }
  }

  Rcpp::CharacterVector res_barcode;
  Rcpp::NumericVector   res_score;

  for(size_t seq1_idx = 0; seq1_idx < reads.size(); seq1_idx++) {
    std::string read = reads[seq1_idx];

    size_t min_distance = UINT_MAX;
    std::string best_barcode;

    for(size_t seq2_idx = 0; seq2_idx < barcodes.size(); seq2_idx++) {
      std::string barcode = barcodes[seq2_idx];

      unsigned int distance = 0;

      for (size_t i = 0; i < n; i++) {
        if (read[i] != barcode[i]) {
          distance += cost_sub_;
        }

      }
      if (distance < min_distance) {
        best_barcode = barcode;
        min_distance = distance;
      }
    }
    res_barcode.push_back(best_barcode);
    res_score.push_back(min_distance);
  }
  return Rcpp::DataFrame::create( Named("barcode")= res_barcode, Named("distance") = res_score, Named("stringsAsFactors") = false);
}

unsigned int HammingDistance::static_min_set_distance(const std::vector<Sequence> &sequences, const size_t n, const unsigned int cost_sub) {
  size_t n_elements = sequences.size();

  unsigned int set_min_distance = UINT_MAX;

  for (size_t row_idx(0); row_idx < n_elements; row_idx++) {
    for (size_t col_idx(row_idx+1); col_idx < n_elements; col_idx++) {

      Sequence sequence1 = sequences.at(row_idx);
      Sequence sequence2 = sequences.at(col_idx);
      
      unsigned int min_distance  = static_distance(sequence1, sequence2, cost_sub);

      if (min_distance < set_min_distance)
        set_min_distance = min_distance;
    }
  }
  return(set_min_distance);
}

unsigned int HammingDistance::static_min_seq_distance(const std::vector< Sequence > &seqs, const Sequence &sequence, const size_t n, const unsigned int cost_sub) {
  size_t n_elements = seqs.size();

  unsigned int global_min_dist = UINT_MAX;
  
  for (size_t seq1_index = 0; seq1_index < n_elements; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance  = static_distance(sequence1,sequence, cost_sub);

    if (min_distance < global_min_dist)
      global_min_dist = min_distance;

  }

  return(global_min_dist);
}

unsigned int HammingDistance::static_distance(const Sequence &sequence1, const Sequence &sequence2, const unsigned int cost_sub) {
  size_t n = sequence1.length();
  size_t m = sequence2.length();

  if (n != m) {
    // Really bad! Bad dog!
    n = std::min(n,m);
  }

  unsigned int min_distance = 0;

  for (size_t i = 0; i < n; i++) {
    if (sequence1.at(i) != sequence2.at(i)) {
      min_distance += cost_sub;
    }
  }

  return (min_distance);
}


