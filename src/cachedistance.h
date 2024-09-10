// cachedistance.h

#ifndef __CACHEDISTANCE_H_INCLUDED__ 
#define __CACHEDISTANCE_H_INCLUDED__

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include "distance.h"

#include <memory>
#include <map>

// Boost symmetric matrix
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/shared_ptr.hpp>

class CacheKey {
  public:

    CacheKey(const Sequence &s1, const Sequence &s2);

    Sequence s1_, s2_;

    bool operator<(const CacheKey& ) const;

};

class CacheDistance : public Distance {
  public:
    CacheDistance(const boost::shared_ptr< Distance > &dist); // no other constructor
    virtual ~CacheDistance();

    virtual unsigned int distance(const Sequence &seq1, const Sequence &seq2);
    virtual unsigned int min_seq_distance(  const std::vector< Sequence > &, const Sequence &, const size_t);
    virtual unsigned int min_set_distance(  const std::vector< Sequence > &, const size_t);
    virtual bool         is_seq_insertable( const std::vector< Sequence > &, const Sequence &, const size_t n, const unsigned int min_dist);
    virtual Rcpp::DataFrame demultiplex(const std::vector< std::string > &, const std::vector< std::string > &);

  private:
    CacheDistance(); // no other constructor
    boost::shared_ptr<Distance> dist_;
    std::map< CacheKey, unsigned int > cache_;

};

#endif 
