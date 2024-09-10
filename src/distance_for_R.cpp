// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "distance.h"
#include "sequencelevenshteindistance.h"
#include "levenshteindistance.h"
#include "hammingdistance.h"
#include "phaseshiftdist.h"

#include "create_distance_func.h"

using namespace Rcpp;

// [[Rcpp::export(".distance")]]
unsigned long int distance(std::string sequence1, std::string sequence2, std::string metric, unsigned int cost_sub, unsigned int cost_indel) {
  // Get the right distance algorithm
  boost::shared_ptr< Distance > dist = create_distance_func(metric, cost_sub, cost_indel);

  return dist->distance(Sequence(sequence1), Sequence(sequence2));
}
