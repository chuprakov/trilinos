// ////////////////////////////////////////////////////////////////////
// qp_solver_stats.h

#ifndef QP_SOLVER_STATS_HH
#define QP_SOLVER_STATS_HH

#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackTypes.h"
#include "ConstrainedOptimizationPack/include/QPSolverStats.h"
#include "GeneralIterationPack/include/CastIQMember.h"

namespace ReducedSpaceSQPPack {

/// Name given to the active set statistics iteration quantity
extern const std::string qp_solver_stats_name;

using ConstrainedOptimizationPack::QPSolverStats;

///
/** Class for object that attempts to return an IterQuantityAccess<QPSolverStats>
  * from an AlgorithmState object with the name qp_solver_stats_name.
  */
class qp_solver_stats_iq_member
	: public CastIQMember<QPSolverStats>
{
public:
    qp_solver_stats_iq_member()
    	: CastIQMember<QPSolverStats>(qp_solver_stats_name)
    {}
};

}	// end namespace ReducedSpaceSQPPack

#endif	// QP_SOLVER_STATS_HH
