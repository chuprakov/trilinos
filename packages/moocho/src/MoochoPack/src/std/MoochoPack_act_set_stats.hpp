// ////////////////////////////////////////////////////////////////////
// act_set_stats.h

#ifndef ACT_SET_STATS_HH
#define ACT_SET_STATS_HH

#include "ActSetStats.h"

namespace ReducedSpaceSQPPack {

/// Name given to the active set statistics iteration quantity
extern const std::string act_set_stats_name;

///
/** Function attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name act_set_stats_name.
  */
IterQuantityAccess<ActSetStats>& act_set_stats( AlgorithmState& state );

///
const IterQuantityAccess<ActSetStats>& act_set_stats( const AlgorithmState& state );


}	// end namespace ReducedSpaceSQPPack

#endif	// ACT_SET_STATS_HH