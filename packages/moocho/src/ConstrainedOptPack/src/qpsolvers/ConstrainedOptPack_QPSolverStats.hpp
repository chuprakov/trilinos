// ////////////////////////////////////////////////////////////////////
// QPSolverStats.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef COP_QP_SOLVER_STATS_H
#define COP_QP_SOLVER_STATS_H

#include "ConstrainedOptPack/src/ConstrainedOptPackTypes.hpp"

namespace ConstrainedOptPack {

///
/** Class for storing statistics about a run of a (active set?) QP solver.
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
	/// Enumeration for the type of projected QP on output
	enum EConvexity {
		CONVEXITY_NOT_KNOWN = static_cast<int>(NOT_KNOWN),
		CONVEX              = 0,
		NONCONVEX           = 1
	};

	// Public interface

	/// Construct all unknowns
	QPSolverStats()
		: solution_type_(SOLUTION_TYPE_NOT_KNOWN)
		, convexity_(CONVEXITY_NOT_KNOWN)
		, num_qp_iter_(NOT_KNOWN)
		, num_adds_(NOT_KNOWN), num_drops_(NOT_KNOWN)
		, warm_start_(false), infeasible_qp_(false)
	{}
	/// Initialize the statistics
	void set_stats(
		ESolutionType solution_type, EConvexity convexity
		,int num_qp_iter, int num_adds, int num_drops
		, bool warm_start, bool infeasible_qp )
	{
		solution_type_	= solution_type;
		convexity_      = convexity;
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
	EConvexity convexity() const
	{
		return convexity_;
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
	EConvexity      convexity_;
	int				num_qp_iter_;
	int				num_adds_;
	int				num_drops_;
	bool			warm_start_;
	bool			infeasible_qp_;

};	// end class QPSolverStats

inline
std::string toString( const QPSolverStats::ESolutionType &solution_type )
{
	switch(solution_type) {
		case QPSolverStats::SOLUTION_TYPE_NOT_KNOWN:
			return "SOLUTION_TYPE_NOT_KNOWN";
			break;
		case QPSolverStats::OPTIMAL_SOLUTION:
			return "OPTIMAL_SOLUTION";
			break;
		case QPSolverStats::PRIMAL_FEASIBLE_POINT:
			return "PRIMAL_FEASIBLE_POINT";
			break;
		case QPSolverStats::DUAL_FEASIBLE_POINT:
			return "DUAL_FEASIBLE_POINT";
			break;
		case QPSolverStats::SUBOPTIMAL_POINT:
			return "SUBOPTIMAL_POINT";
			break;
		default:
			TEST_FOR_EXCEPT(true);
	}
}

}	// end namespace ConstrainedOptPack

#endif	// COP_QP_SOLVER_STATS_H
