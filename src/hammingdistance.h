// hammingdistance.h

#ifndef __HAMMINGDISTANCE_H_INCLUDED__ 
#define __HAMMINGDISTANCE_H_INCLUDED__

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "distance.h"

#include <memory>

#define HAMMINGDISTANCE_DEFAULT_COST_SUB 1
class HammingDistance : public Distance {

  public:
    HammingDistance(const unsigned int cost_sub = HAMMINGDISTANCE_DEFAULT_COST_SUB); // default constructor
    virtual ~HammingDistance();

    virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2);
    virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t);
    virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t);
    virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist);
    virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &);
    
    static unsigned int static_distance(const Sequence &, const Sequence &, const unsigned int cost_sub = HAMMINGDISTANCE_DEFAULT_COST_SUB);
    static unsigned int static_min_seq_distance(const std::vector< Sequence > &, const Sequence &, const size_t, const unsigned int cost_sub = HAMMINGDISTANCE_DEFAULT_COST_SUB);
    static unsigned int static_min_set_distance(const std::vector< Sequence > &, const size_t, const unsigned int cost_sub = HAMMINGDISTANCE_DEFAULT_COST_SUB);

  private:
    unsigned int cost_sub_;

};

#endif 
