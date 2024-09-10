// in distance.h

#ifndef __DISTANCE_H_INCLUDED__ 
#define __DISTANCE_H_INCLUDED__

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include <vector>
#include <memory>
#include <string>

#include "sequence.h"

class Distance {
   public:
      // pure virtual function
      virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2) = 0;
      virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t) = 0;
      virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t) = 0;
      virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist) = 0;
      virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &) = 0;
};

#endif

