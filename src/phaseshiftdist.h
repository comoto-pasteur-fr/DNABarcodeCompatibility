// phaseshiftdist.h

#ifndef __PHASESHIFTDIST_H_INCLUDED__ 
#define __PHASESHIFTDIST_H_INCLUDED__

#include "distance.h"

#include <memory>

#include <boost/shared_ptr.hpp>

#define PHASESHIFT_DEFAULT_COST_SHIFT 1
#define PHASESHIFT_DEFAULT_COST_SUB   1

typedef std::pair<Sequence, unsigned int > dist_pair_t;

class PhaseshiftDist : public Distance {

  public:
    PhaseshiftDist(const unsigned int cost_sub = PHASESHIFT_DEFAULT_COST_SUB, const unsigned int cost_shift = PHASESHIFT_DEFAULT_COST_SHIFT); 
    virtual ~PhaseshiftDist();

    virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2);
    virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t);
    virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t);
    virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist);
    virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &);
    
    static unsigned int static_distance(const Sequence&, const Sequence&, const unsigned int cost_sub = PHASESHIFT_DEFAULT_COST_SUB, const unsigned int cost_shift = PHASESHIFT_DEFAULT_COST_SHIFT);
    static unsigned int static_limited_distance(const Sequence &, const Sequence &, const unsigned int max_dist, const unsigned int cost_sub = PHASESHIFT_DEFAULT_COST_SUB, const unsigned int cost_shift = PHASESHIFT_DEFAULT_COST_SHIFT);
    static unsigned int static_min_seq_distance(const std::vector<Sequence> &, const Sequence&, const size_t, const unsigned int cost_sub = PHASESHIFT_DEFAULT_COST_SUB, const unsigned int cost_shift = PHASESHIFT_DEFAULT_COST_SHIFT);
    static unsigned int static_min_set_distance(const std::vector<Sequence> &, const size_t , const unsigned int cost_sub = PHASESHIFT_DEFAULT_COST_SUB, const unsigned int cost_shift = PHASESHIFT_DEFAULT_COST_SHIFT);

  private:
    unsigned int cost_sub_, cost_shift_;
};

struct CompareDistPair
{
      bool operator()(std::pair<Sequence , unsigned int > const& a, std::pair<Sequence , unsigned int > const& b) const { return b.second < a.second; }
};

#endif 
