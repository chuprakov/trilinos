// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinSimpleNewtonSolverDecl.hpp

#ifndef TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_DECL_HPP
#define TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_DECL_HPP

#include "TSFCoreNonlinTypes.hpp"
#include "StandardMemberCompositionMacros.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** A simple newton-based nonlinear equation solver with backtracking
 * line search.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SimpleNewtonSolver {
public:

	///
	/** Set the tolerance on ||c|| <= tol for convergence.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, tol );

	///
	/** Set the maximum number of newton iterations to take.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxNewtonIter );

	///
	/** Set the maximum number of backtracking line search iterations to take.
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( int, maxLineSearchIter );

	///
	/** Set the armijo constant for the line search
	 */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( Scalar, eta );

	///
	SimpleNewtonSolver(
		const Scalar   tol                    = 1e-8 
		,const int     maxNewtonIter          = 1000
		,const int     maxLineSearchIter      = 20
		,const Scalar  eta                    = 1e-4
		);
		
	///
	/** Solve a set of nonlinear equations given an initial guess.
	 *
	 * @param  np  [in/out] Defines the nonlinear problem to be solved.  Note that
	 *             <tt>np->initialize()</tt> must be called prior to calling this function.
	 * @param  y   [in/out] On input, <tt>y</tt> contains the initial guess.  On output, <tt>y</tt>
	 *             contains the solution (or partial solution) of the set of nonlinear equations.
	 *
	 * ToDo: Finish documentation!
	 */
	Solvers::SolveReturn solve(
		NonlinearProblemFirstOrder<Scalar>          *np
		,Vector<Scalar>                             *y
		,std::ostream                               *out     = NULL
		,bool                                       dumpAll  = false
		) const;

}; // class SimpleNewtonSolver

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_DECL_HPP
