// ////////////////////////////////////////////////////////////////////
// qp_solver_stats.h

#ifndef QP_SOLVER_STATS_HH
#define QP_SOLVER_STATS_HH

#include "QPSolverStats.h"

namespace ReducedSpaceSQPPack {

/// Name given to the active set statistics iteration quantity
extern const std::string qp_solver_stats_name;

///
/** Function attempts to return an IterQuantityAccess<ActSetStats>
  * from an AlgorithmState object with the name qp_solver_stats_name.
  */
IterQuantityAccess<QPSolverStats>& qp_solver_stats( AlgorithmState& state );

///
const IterQuantityAccess<QPSolverStats>& qp_solver_stats( const AlgorithmState& state );

}	// end namespace ReducedSpaceSQPPack

#endif	// QP_SOLVER_STATS_HH