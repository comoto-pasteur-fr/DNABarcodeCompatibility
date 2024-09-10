// in phaseshiftdist.cpp
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sequence.h"
#include "phaseshiftdist.h"

#include <assert.h>
#include <limits.h>

#include <algorithm>
#include <memory>
#include <queue>
#include <utility>
#include <set>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/cstdint.hpp>

using namespace Rcpp;

PhaseshiftDist::PhaseshiftDist(const unsigned int cost_sub, const unsigned int cost_shift) : cost_sub_(cost_sub), cost_shift_(cost_shift) {

}

PhaseshiftDist::~PhaseshiftDist() {};

unsigned int PhaseshiftDist::distance(const Sequence &sequence1, const Sequence &sequence2) {
  return(static_limited_distance(sequence1, sequence2, std::numeric_limits<unsigned int>::max(), cost_sub_, cost_shift_));
}


    
unsigned int  PhaseshiftDist::min_seq_distance(const std::vector<Sequence> &seqs, const Sequence& sequence, const size_t n) {
  return(static_min_seq_distance( seqs, sequence, n, cost_sub_, cost_shift_));
}


unsigned int  PhaseshiftDist::min_set_distance(const std::vector<Sequence> &sequences, const size_t n) {
  return(static_min_set_distance(sequences, n, cost_sub_, cost_shift_));
}


bool PhaseshiftDist::is_seq_insertable(const std::vector<Sequence> &seqs, const Sequence &sequence, const size_t n, const unsigned int min_dist) {
  size_t n_elements = seqs.size();
  bool seq_is_insertable = true;

  for (size_t seq1_index = 0; seq1_index < n_elements && seq_is_insertable; seq1_index++) {
    //std::cerr << ".";
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance = static_limited_distance(sequence1, sequence, min_dist, cost_sub_, cost_shift_);

    if (min_distance < min_dist)
      seq_is_insertable = false;

  }

  return seq_is_insertable;
}

Rcpp::DataFrame PhaseshiftDist::demultiplex(const std::vector< std::string > &barcodes, const std::vector< std::string > &reads) {
  if ((barcodes.size() < 2) || (reads.size() < 1)) {
    Rcpp::stop("At least one read and two barcodes need to be provided");
  }

  size_t n = barcodes[1].length();

  for (size_t i = 1; i < barcodes.size(); i++) {
    if (barcodes[i].length() != n) {
      Rcpp::stop("Length of all barcodes must be equal. (At the moment).");
    }
  }

  Rcpp::CharacterVector res_barcode;
  Rcpp::NumericVector   res_score;

  for(size_t seq1_idx = 0; seq1_idx < reads.size(); seq1_idx++) {
    Sequence read_seq(reads[seq1_idx]);

    std::string best_barcode;
    unsigned int min_distance = UINT_MAX;

    for(size_t seq2_idx = 0; seq2_idx < barcodes.size(); seq2_idx++) {
      std::string barcode = barcodes[seq2_idx];
      Sequence barcode_seq(barcode);

      unsigned int distance = this->distance(read_seq, barcode_seq);

      if (distance < min_distance) {
        best_barcode = barcode;
        min_distance = distance;
      }
    }
    //Rprintf("Best: %s Dist: %i\n", best_barcode.c_str(), min_distance);
    res_barcode.push_back(best_barcode);
    res_score.push_back(min_distance);
  }

  return Rcpp::DataFrame::create( Named("barcode")= res_barcode, Named("distance") = res_score, Named("stringsAsFactors") = false);
}

unsigned int PhaseshiftDist::static_distance(const Sequence &sequence1, const Sequence& sequence2, const unsigned int cost_sub, const unsigned int cost_shift) {
  return(static_limited_distance(sequence1, sequence2, std::numeric_limits<unsigned int>::max(), cost_sub, cost_shift));
}

unsigned int PhaseshiftDist::static_limited_distance(const Sequence &sequence1, const Sequence &sequence2, const unsigned int max_dist, const unsigned int cost_sub, const unsigned int cost_shift) {

  if (sequence1 == sequence2)
    return 0;

  // Get local copies to circumvent shared_ptr problems

  size_t n = sequence1.length();
  size_t m = sequence2.length();

  if (n != m) Rcpp::stop("Lengths of both sequences need to be equal.");

  // Woah
  std::priority_queue<dist_pair_t, std::vector<dist_pair_t>, CompareDistPair > dist_queue;

  // FIXME: replace by unordered set
  std::set<Sequence> already_seen;

  dist_queue.push(dist_pair_t(sequence1, 0));
  already_seen.insert(sequence1);

  unsigned int curr_dist = std::min(n * cost_sub, n * cost_shift);

  while (!dist_queue.empty()) {
    // Retrieve head of queue
    dist_pair_t dist_pair = dist_queue.top(); dist_queue.pop();

    // Extract sequence of head
    Sequence mutated_seq    = dist_pair.first;
    // Extract current distance of head from the compared sequence
    unsigned int mutations  = dist_pair.second;

    //Rcpp::Rcerr << "Testing " << *mutated_seq << " (d=" << mutations << ")" << std::endl;

    /***********************
    * Test return criteria *
    ************************/

    // Cut short, if this seq is already too far away from seq2
    if (mutations >= max_dist)
      return max_dist;

    // Return, if we already had a sequence that is closer to seq2 (through phase shifts + substitutions)
    if (mutations >= curr_dist)
      return curr_dist;

    // Check if we this sequence and sequence 2 are equal (end of function)
    if (mutated_seq == sequence2)
      return mutations;
                
    /************************************************************
    * Calculate Hamming distance after having finished shifting *
    *************************************************************/

    unsigned int h_dist = 0;

    // only go to the end of the shorter sequence
    size_t min_n = std::min(mutated_seq.length(),sequence2.length());
    for (size_t i = 0; i < min_n; i++) {
      if (mutated_seq.at(i) != sequence2.at(i)) {
        h_dist += cost_sub;
      }
    }

    //Rcpp::Rcerr << "Total is "  << mutations + h_dist << std::endl;

    /*********************
    * Find new benchmark *
    **********************/
    // Set a new benchmark for smallest distance
    if (mutations + h_dist < curr_dist) {
      curr_dist = mutations + h_dist;
    }

    /***************************************
    * Try frontal insertions and deletions *
    ****************************************/

    if (!(mutations + cost_shift > curr_dist) and !(mutations + cost_shift > max_dist)) {
      // temporary sequence for new mutations
      //Sequence tmp_seq;

      // Try frontal insertions
      BOOST_FOREACH(uint64_t ins_base, Sequence::REAL_BASES) {
        // Insertion
        Sequence tmp_seq = mutated_seq.insert(0, ins_base);

        if (already_seen.find(tmp_seq) == already_seen.end()) {
          already_seen.insert(tmp_seq);
          dist_queue.push(dist_pair_t(tmp_seq, cost_shift+mutations));
        }

      }
      // Try frontal deletion
      Sequence tmp_seq = mutated_seq.remove(0);
      if (already_seen.find(tmp_seq) == already_seen.end()) {
        already_seen.insert(tmp_seq);
        dist_queue.push(dist_pair_t(tmp_seq, cost_shift+mutations));
      }
    }
  }

  return curr_dist; // Should never every be called
}

unsigned int  PhaseshiftDist::static_min_seq_distance(const std::vector<Sequence> &seqs, const Sequence &sequence, const size_t n, const unsigned int cost_sub, const unsigned int cost_shift) {

  unsigned int global_min_dist = UINT_MAX;
  
  size_t n_elements = seqs.size();
  for (size_t seq1_index = 0; seq1_index < n_elements; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    unsigned int min_distance = static_distance(sequence1, sequence, cost_sub, cost_shift);

    if (min_distance < global_min_dist)
      global_min_dist = min_distance;

  }

  return(global_min_dist);
}

unsigned int  PhaseshiftDist::static_min_set_distance(const std::vector<Sequence> &sequences, const size_t n, const unsigned int cost_sub, const unsigned int cost_shift) {

  unsigned int set_min_distance = UINT_MAX;

  size_t n_elements = sequences.size();
  for (size_t row_idx(0); row_idx < n_elements; row_idx++) {
    for (size_t col_idx(row_idx+1); col_idx < n_elements; col_idx++) {
      Sequence sequence1 = sequences.at(row_idx);
      Sequence sequence2 = sequences.at(col_idx);
      
      unsigned int min_distance = static_distance(sequence1, sequence2, cost_sub, cost_shift);

      if (min_distance < set_min_distance)
        set_min_distance = min_distance;
    }
  }
  return(set_min_distance);
}


