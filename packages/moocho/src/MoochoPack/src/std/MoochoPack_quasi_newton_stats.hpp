// ////////////////////////////////////////////////////////////////////
// quasi_newton_stats.h

#ifndef QUASI_NEWTON_STATS_HH
#define QUASI_NEWTON_STATS_HH

#include "QuasiNewtonStats.h"

namespace ReducedSpaceSQPPack {

/// Name given to the quasi-Newton updating staistics iteration quantity
extern const std::string quasi_newton_stats_name;

///
/** Function attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name quasi_newton_stats_name.
  */
IterQuantityAccess<QuasiNewtonStats>& quasi_newton_stats( AlgorithmState& state );

///
const IterQuantityAccess<QuasiNewtonStats>& quasi_newton_stats( const AlgorithmState& state );


}	// end namespace ReducedSpaceSQPPack

#endif	// QUASI_NEWTON_STATS_HH