// ////////////////////////////////////////////////////////////////////
// qp_solver_stats.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include "../../include/std/qp_solver_stats.h"
#include "GeneralIterationPack/include/cast_iq.h"

extern const std::string ReducedSpaceSQPPack::qp_solver_stats_name = "qp_solver_stats";

// Warning: A MS VC++ compiler error does not allow the proper
// definition of this function using the namespace prefix so we
// have to use the unsafe namespace enclosure.  There is therefore
// no way for the compiler to determine if this definition actually
// defines a function already declared or if it is a new function
// all together (which is the case in the global namespace).

namespace ReducedSpaceSQPPack {

IterQuantityAccess<QPSolverStats>& qp_solver_stats( AlgorithmState& state )
{
	return GeneralIterationPack::cast_iq<QPSolverStats>( state, qp_solver_stats_name );
}

const IterQuantityAccess<QPSolverStats>& qp_solver_stats( const AlgorithmState& state )
{
	return GeneralIterationPack::cast_iq<QPSolverStats>( state, qp_solver_stats_name );
}

}	// end namespace ReducedSpaceSQPPack 