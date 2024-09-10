// sequencelevenshteindistance.h

#ifndef __SEQUENCELEVENSHTEINDISTANCE_H_INCLUDED__ 
#define __SEQUENCELEVENSHTEINDISTANCE_H_INCLUDED__

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include "distance.h"
#include <memory>

#define SEQUENCELEVENSHTEIN_DEFAULT_COST_SUB 1
#define SEQUENCELEVENSHTEIN_DEFAULT_COST_INDEL 1

class SequenceLevenshteinDistance : public Distance {

  public:
    SequenceLevenshteinDistance(const unsigned int cost_sub = SEQUENCELEVENSHTEIN_DEFAULT_COST_SUB, const unsigned int cost_indel = SEQUENCELEVENSHTEIN_DEFAULT_COST_INDEL);
    virtual ~SequenceLevenshteinDistance();

    virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2);
    virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t);
    virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t);
    virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist);
    virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &);
    
    static unsigned int static_distance(const Sequence&, const Sequence&, const unsigned int cost_sub = SEQUENCELEVENSHTEIN_DEFAULT_COST_SUB, const unsigned int cost_indel = SEQUENCELEVENSHTEIN_DEFAULT_COST_INDEL);
    static unsigned int static_min_seq_distance(const std::vector<Sequence> &, const Sequence &, const size_t, const unsigned int cost_sub = SEQUENCELEVENSHTEIN_DEFAULT_COST_SUB, const unsigned int cost_indel = SEQUENCELEVENSHTEIN_DEFAULT_COST_INDEL);
    static unsigned int static_min_set_distance(const std::vector<Sequence> &, const size_t, const unsigned int cost_sub = SEQUENCELEVENSHTEIN_DEFAULT_COST_SUB, const unsigned int cost_indel = SEQUENCELEVENSHTEIN_DEFAULT_COST_INDEL);

  private:
    unsigned int cost_sub_, cost_indel_;

};

#endif 
