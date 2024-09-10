// levenshteindistance.h

#ifndef __LEVENSHTEINDISTANCE_H_INCLUDED__ 
#define __LEVENSHTEINDISTANCE_H_INCLUDED__

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "distance.h"

#include <memory>

#define LEVENSHTEINDISTANCE_DEFAULT_COST_SUB 1
#define LEVENSHTEINDISTANCE_DEFAULT_COST_INDEL 1

class LevenshteinDistance : public Distance {

  public:
    LevenshteinDistance(const unsigned int cost_sub = LEVENSHTEINDISTANCE_DEFAULT_COST_SUB, const unsigned int cost_indel = LEVENSHTEINDISTANCE_DEFAULT_COST_INDEL); 
    virtual ~LevenshteinDistance();

    virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2);
    virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t);
    virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t);
    virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist);
    virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &);

    static unsigned int static_distance(const Sequence&, const Sequence&, const unsigned int, const unsigned int);
    static unsigned int static_min_seq_distance(const std::vector<Sequence> &, const Sequence &, const size_t, const unsigned int, const unsigned int);
    static unsigned int static_min_set_distance(const std::vector<Sequence> &, const size_t, const unsigned int, const unsigned int);

  private:
    unsigned int cost_sub_, cost_indel_;
};

#endif 
