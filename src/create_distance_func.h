#include "distance.h"

#include <boost/shared_ptr.hpp>

boost::shared_ptr< Distance > create_distance_func(const std::string metric, const unsigned int cost_sub, const unsigned int cost_indel);
