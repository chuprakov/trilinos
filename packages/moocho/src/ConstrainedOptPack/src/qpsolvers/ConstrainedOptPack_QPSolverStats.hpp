// ////////////////////////////////////////////////////////////////////
// QPSolverStats.h

#ifndef COP_QP_SOLVER_STATS_H
#define COP_QP_SOLVER_STATS_H

#include "ConstrainedOptimizationPackTypes.h"

namespace ConstrainedOptimizationPack {

///
/** Class for storing statistics about the active-set QP solver.
  */
class QPSolverStats {
public:

	// Public types

	/// Set to this value if a statistic is not known.
	enum { NOT_KNOWN = -1 };

	/// Enumeration for the type of point returned from solve_qp(...).
	enum ESolutionType {
		SOLUTION_TYPE_NOT_KNOWN = static_cast<int>(NOT_KNOWN),
		OPTIMAL_SOLUTION		= 0,
		PRIMAL_FEASIBLE_POINT	= 1,
		DUAL_FEASIBLE_POINT		= 2,
		SUBOPTIMAL_POINT		= 3
		};

	// Public interface

	/// Construct all unknowns
	QPSolverStats()
		: solution_type_(SOLUTION_TYPE_NOT_KNOWN), num_qp_iter_(NOT_KNOWN)
			, num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
			, warm_start_(false), infeasible_qp_(false)
	{}

	/// Initialize the statistics
	void set_stats( ESolutionType solution_type, int num_qp_iter
		, int num_adds, int num_drops, bool warm_start, bool	infeasible_qp )
	{
		solution_type_	= solution_type;
		num_qp_iter_	= num_qp_iter; 
		num_adds_		= num_adds;
		num_drops_		= num_drops;
		warm_start_		= warm_start;
		infeasible_qp_	= infeasible_qp;
	}

	///
	ESolutionType solution_type() const
	{
		return solution_type_;
	}

	///
	int num_qp_iter() const
	{
		return num_qp_iter_;
	}
	///
	int	num_adds() const
	{
		return num_adds_;
	}
	///
	int	num_drop() const
	{
		return num_drops_;
	}
	///
	int	warm_start() const
	{
		return warm_start_;
	}
	///
	int	infeasible_qp() const
	{
		return infeasible_qp_;
	}

private:
	ESolutionType	solution_type_;
	int				num_qp_iter_;
	int				num_adds_;
	int				num_drops_;
	bool			warm_start_;
	bool			infeasible_qp_;

};	// end class ActSetStats

}	// end namespace ConstrainedOptimizationPack

#endif	// COP_QP_SOLVER_STATS_H