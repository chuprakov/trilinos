// ////////////////////////////////////////////////////////////////////
// act_set_stats.h

#ifndef ACT_SET_STATS_HH
#define ACT_SET_STATS_HH

#include "ActSetStats.h"
#include "GeneralIterationPack/include/CastIQMember.h"

namespace ReducedSpaceSQPPack {

/// Name given to the active set statistics iteration quantity
extern const std::string act_set_stats_name;

///
/** Class for object that attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name act_set_stats_name.
  */
class act_set_stats_iq_member
	: public CastIQMember<ActSetStats>
{
public:
    act_set_stats_iq_member()
    	: CastIQMember<ActSetStats>(act_set_stats_name)
    {}
};

}	// end namespace ReducedSpaceSQPPack

#endif	// ACT_SET_STATS_HH