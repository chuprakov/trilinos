// ///////////////////////////////////////////////////////////////////////
// ReducedSpaceSQPPackExceptions.h

#ifndef REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
#define REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H

#include "ReducedSpaceSQPPackTypes.h"
#include "ConstrainedOptimizationPack/include/QPSolverStats.h"

namespace ReducedSpaceSQPPack {

// Thrown if the constraints are infeasible
class InfeasibleConstraints : public std::logic_error
{public: InfeasibleConstraints(const std::string& what_arg) : std::logic_error(what_arg) {}};

/// Thrown if a line search failure occurs.
class LineSearchFailure : public std::runtime_error
{public: LineSearchFailure(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a runtime test failed.
class TestFailed : public std::runtime_error
{public: TestFailed(const std::string& what_arg) : std::runtime_error(what_arg){}};

/// Thrown if a the QP failed and was not corredted
class QPFailure : public std::runtime_error
{
public:
	QPFailure(const std::string& what_arg
			  , const ConstrainedOptimizationPack::QPSolverStats& _qp_stats)
		: std::runtime_error(what_arg)
		, qp_stats(_qp_stats)
		{}
	ConstrainedOptimizationPack::QPSolverStats qp_stats;
};

}	// end namespace ReducedSpaceSQPPack 

#endif // REDUCED_SPACE_SQP_PACK_EXCEPTIONS_H
