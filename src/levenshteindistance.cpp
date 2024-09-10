// in sequencelevenshteindistance.cpp

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "sequence.h"
#include "levenshteindistance.h"

#include <limits.h>

#include <algorithm>
#include <memory>

using namespace Rcpp;

LevenshteinDistance::LevenshteinDistance(const unsigned int cost_sub, const unsigned int cost_indel) : cost_sub_(cost_sub), cost_indel_(cost_indel)  {

}

LevenshteinDistance::~LevenshteinDistance() {};

unsigned int LevenshteinDistance::distance(const Sequence &seq1, const Sequence &seq2) {
    return(static_distance(seq1, seq2, cost_sub_, cost_indel_));
}

unsigned int LevenshteinDistance::min_seq_distance(const std::vector< Sequence > &seqs, const Sequence &sequence, const size_t n) {
    return(static_min_seq_distance(seqs, sequence, n, cost_sub_, cost_indel_));
}

unsigned int LevenshteinDistance::min_set_distance(const std::vector<Sequence> &sequences, const size_t n) {
    return(static_min_set_distance(sequences, n, cost_sub_, cost_indel_));
}


bool LevenshteinDistance::is_seq_insertable(const std::vector<Sequence> &seqs, const Sequence &sequence, const size_t n, const unsigned int min_dist) {
  size_t n_elements = seqs.size();

  size_t m = sequence.length();

  unsigned int distances[n+1][m+1];

  for (size_t i = 0; i <= n; i++) {
    distances[i][0] = i * cost_indel_;
  }
  
  for (size_t j = 0; j <= m; j++) {
    distances[0][j] = j * cost_indel_;
  }

  bool seq_is_insertable = true;

  for (size_t seq1_index = 0; seq1_index < n_elements && seq_is_insertable; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    for (size_t i = 1; i <= n; i++) {
      for (size_t j = 1; j <= m; j++) {
        unsigned int cost = 0;

        if (sequence1.at(i-1) != sequence.at(j-1)) {
          cost = cost_sub_;
        }

        distances[i][j] = std::min(std::min(distances[i-1][j-1] + cost, distances[i][j-1] + cost_indel_), distances[i-1][j] + cost_indel_);

      }
    }

    unsigned int min_distance  = distances[n][m];

    if (min_distance < min_dist)
      seq_is_insertable = false;
  }

  return seq_is_insertable;

}

Rcpp::DataFrame LevenshteinDistance::demultiplex(const std::vector< std::string > &barcodes, const std::vector< std::string > &reads) {
  if ((barcodes.size() < 2) || (reads.size() < 1)) {
    Rcpp::stop("At least one read and two barcodes need to be provided");
  }

  size_t n = barcodes[1].length();

  for (unsigned int i = 1; i < barcodes.size(); i++) {
    if (barcodes[i].length() != n) {
      Rcpp::stop("Length of all barcodes must be equal. (At the moment).");
    }
  }

  Rcpp::CharacterVector res_barcode;
  Rcpp::NumericVector   res_score;

  for(unsigned int read_idx = 0; read_idx < reads.size(); read_idx++) {
    // Extract read of indexes 0..j..m-1
    std::string read = reads[read_idx];
    // Get length of read
    size_t m = read.length();

    // Prepare an array to contain the distances of the dynamic algorithm
    unsigned int distances[n+1][m+1];
    // Prepare distances of dyn algo with distance of empty string to positions
    for (size_t i = 0; i <= n; i++) {
      distances[i][0] = i * cost_indel_;
    }
    for (size_t j = 0; j <= m; j++) {
      distances[0][j] = j * cost_indel_;
    }

    // Best distance yet  
    unsigned int min_distance = UINT_MAX;
    // Best barcode yet
    std::string best_barcode;

    for(unsigned int barcode_idx = 0; barcode_idx < barcodes.size(); barcode_idx++) {
      // Extract barcode of indexes 0..i..n-1
      std::string barcode = barcodes[barcode_idx];

      // Dyn Algo
      for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= m; j++) {
          unsigned int cost = 0;

          if (barcode[i-1] != read[j-1]) {
            cost = cost_sub_;
          }

          distances[i][j] = std::min(std::min(distances[i-1][j-1] + cost, distances[i][j-1] + cost_indel_), distances[i-1][j] + cost_indel_);

        }
      }

      unsigned int distance = distances[n][m];

      // Update best found barcode?
      if (distance < min_distance) {
        min_distance = distance;
        best_barcode = barcode;
      }
    }
    res_barcode.push_back(best_barcode);
    res_score.push_back(min_distance);
  }
  return Rcpp::DataFrame::create( Named("barcode")= res_barcode, Named("distance") = res_score, Named("stringsAsFactors") = false);
}

unsigned int LevenshteinDistance::static_min_set_distance(const std::vector<Sequence> &sequences, const size_t n, const unsigned int cost_sub, const unsigned int cost_indel) {
  size_t n_elements = sequences.size();

  //unsigned int multi_distances[n_elements][n_elements];

  unsigned int set_min_distance = UINT_MAX;

  unsigned int distances[n+1][n+1];

  for (size_t i = 0; i <= n; i++) {
    distances[i][0] = i * cost_indel;
  }
  for (size_t j = 0; j <= n; j++) {
    distances[0][j] = j * cost_indel;
  }
  
  for (size_t row_idx(0); row_idx < n_elements; row_idx++) {
    for (size_t col_idx(row_idx+1); col_idx < n_elements; col_idx++) {

      Sequence sequence1 = sequences.at(row_idx);
      Sequence sequence2 = sequences.at(col_idx);
      
      for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= n; j++) {
          unsigned int cost = 0;

          if (sequence1.at(i-1) != sequence2.at(j-1)) {
            cost = cost_sub;
          }

          distances[i][j] = std::min(std::min(distances[i-1][j-1] + cost, distances[i][j-1] + cost_indel), distances[i-1][j] + cost_indel);

        }
      }

      unsigned int min_distance  = distances[n][n];

      if (min_distance < set_min_distance)
        set_min_distance = min_distance;
    }
  }
  return(set_min_distance);
}

unsigned int LevenshteinDistance::static_min_seq_distance(const std::vector< Sequence > &seqs, const Sequence &sequence, const size_t n, const unsigned int cost_sub, const unsigned int cost_indel) {

  size_t n_elements = seqs.size();

  size_t m = sequence.length();

  unsigned int global_min_dist = UINT_MAX;
  
  unsigned int distances[n+1][m+1];

  for (size_t i = 0; i <= n; i++) {
    distances[i][0] = i * cost_indel;
  }
  for (size_t j = 0; j <= m; j++) {
    distances[0][j] = j * cost_indel;
  }

  for (size_t seq1_index = 0; seq1_index < n_elements; seq1_index++) {
    Sequence sequence1 = seqs.at(seq1_index);

    for (size_t i = 1; i <= n; i++) {
      for (size_t j = 1; j <= m; j++) {
        unsigned int cost = cost_sub;

        if (sequence1.at(i-1) != sequence.at(j-1)) {
          cost = cost_sub;
        }

        distances[i][j] = std::min(std::min(distances[i-1][j-1] + cost, distances[i][j-1] + cost_indel), distances[i-1][j] + cost_indel);

      }
    }

    unsigned int min_distance  = distances[n][m];

    if (min_distance < global_min_dist)
      global_min_dist = min_distance;

  }

  return(global_min_dist);
}

unsigned int LevenshteinDistance::static_distance(const Sequence &sequence1, const Sequence &sequence2, const unsigned int cost_sub, const unsigned int cost_indel) {
  size_t n = sequence1.length();
  size_t m = sequence2.length();

  unsigned int distances[n+1][m+1];

  // Distance from empty string to the position
  for (size_t i = 0; i <= n; i++) {
    distances[i][0] = i * cost_indel;
  }
  for (size_t j = 0; j <= m; j++) {
    distances[0][j] = j * cost_indel;
  }

  for (size_t i = 1; i <= n; i++) {
    for (size_t j = 1; j <= m; j++) {
      unsigned int cost = 0;

      if (sequence1.at(i-1) != sequence2.at(j-1)) {
        cost = cost_sub;
      }

      distances[i][j] = std::min(std::min(distances[i-1][j-1] + cost, distances[i][j-1] + cost_indel), distances[i-1][j] + cost_indel);
    }
  }

  unsigned int min_distance = distances[n][m];

  return (min_distance);

}


