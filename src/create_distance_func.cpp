// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include "distance.h"
#include "sequencelevenshteindistance.h"
#include "levenshteindistance.h"
#include "hammingdistance.h"
#include "cachedistance.h"
#include "phaseshiftdist.h"

boost::shared_ptr< Distance > create_distance_func(const std::string metric, const unsigned int cost_sub, const unsigned int cost_indel) {
  boost::shared_ptr< Distance > tmp_dist;

  if (0 == metric.compare("hamming")) {
    tmp_dist = boost::shared_ptr<Distance>(new HammingDistance(cost_sub));
  } else if (0 == metric.compare("seqlev")) {
    tmp_dist = boost::shared_ptr<Distance>(new SequenceLevenshteinDistance(cost_sub, cost_indel));
  } else if (0 == metric.compare("levenshtein")) {
    tmp_dist = boost::shared_ptr<Distance>(new LevenshteinDistance(cost_sub, cost_indel));
  } else if (0 == metric.compare("phaseshift")) {
    tmp_dist = boost::shared_ptr<Distance>(new PhaseshiftDist(cost_sub, cost_indel));
  } else {
    Rcpp::stop("Unrecognized distance metric given.");
  }

  return tmp_dist;
}


