// ////////////////////////////////////////////////////////////////////
// quasi_newton_stats.h

#ifndef QUASI_NEWTON_STATS_HH
#define QUASI_NEWTON_STATS_HH

#include "QuasiNewtonStats.h"
#include "GeneralIterationPack/include/CastIQMember.h"

namespace ReducedSpaceSQPPack {

/// Name given to the quasi-Newton updating staistics iteration quantity
extern const std::string quasi_newton_stats_name;

///
/** Class for object that attempts to return an IterQuantityAccess<QuasiNewtonStats>
  * from an AlgorithmState object with the name quasi_newton_stats_name.
  */
class quasi_newton_stats_iq_member
	: public CastIQMember<QuasiNewtonStats>
{
public:
    quasi_newton_stats_iq_member()
    	: CastIQMember<QuasiNewtonStats>(quasi_newton_stats_name)
    {}
};

}	// end namespace ReducedSpaceSQPPack

#endif	// QUASI_NEWTON_STATS_HH