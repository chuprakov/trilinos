// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinSimpleNewtonSolver.hpp

#ifndef TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP
#define TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP

#include "TSFCoreNonlinSimpleNewtonSolverDecl.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreNonlinNonlinearProblem.hpp"
#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "AbstractFactory.hpp"

namespace TSFCore {
namespace Nonlin {

template<class Scalar>
SimpleNewtonSolver<Scalar>::SimpleNewtonSolver(
	const Scalar   tol
	,const int     maxNewtonIter
	,const int     maxLineSearchIter
	,const Scalar  eta
	)
	:tol_(tol), maxNewtonIter_(maxNewtonIter), maxLineSearchIter_(maxLineSearchIter), eta_(eta)
{}
		
template<class Scalar>
Solvers::SolveReturn
SimpleNewtonSolver<Scalar>::solve( NonlinearProblemFirstOrder<Scalar> *np
								   ,Vector<Scalar> *y_inout, std::ostream *out, bool dumpAll ) const
{
	namespace mmp = MemMngPack;
	if(out) *out << "\n*** Entering SimpleNewtonSolver::solve(...) ...\n";
	// Create the data for the problem
	const Index                                     n        = np->space_c()->dim();
	mmp::ref_count_ptr<Vector<Scalar> >             y        = mmp::rcp(y_inout,false);       // Initial guess/solution
	mmp::ref_count_ptr<Vector<Scalar> >             c        = np->space_c()->createMember(); // Residual
	mmp::ref_count_ptr<LinearOpWithSolve<Scalar> >  DcDy     = np->factory_DcDy()->create();  // Jacobian object
	ETransp                                         opDcDy   = np->opDcDy();                  // Transpose argument for DcDy
	mmp::ref_count_ptr<Vector<Scalar> >             dy       = np->space_y()->createMember(); // Newton step for y
	mmp::ref_count_ptr<Vector<Scalar> >             y_new    = np->space_y()->createMember(); // Trial point for y
	// Compute the initial starting point
	np->unsetQuantities(); np->set_c(c.get()); np->set_DcDy(DcDy.get()); // These pointers will be maintained throughout
	np->calc_DcDy(*y,NULL); np->calc_c(*y,NULL,false);
	// Print the starting point
	if(out && dumpAll) {
		*out << "\nInitial starting point:\n";
		*out << "\ny    =\n" << *y   ;
		*out << "\nc    =\n" << *c   ;
		*out << "\nDcDy =\n" << *DcDy;
	}
	// Peform the Newton iterations
	int newtonIter;
	for( newtonIter = 1; newtonIter <= maxNewtonIter(); ++newtonIter ) {
		if(out) *out << "\n*** newtonIter = " << newtonIter << std::endl;
		// Check convergence
		if(out) *out << "\nChecking for convergence ... : ";
		const Scalar phi = 0.5*dot(*c,*c)/n, sqrt_phi = sqrt(phi); // merit function: phi(c) = 0.5*dot(c,c)/n
		const bool isConverged = sqrt_phi <= tol();
		if(out) *out << "sqrt(phi) = sqrt(0.5*dot(c,c)/n) = " << sqrt_phi << ( isConverged ? " <= " : " > " ) << "tol = " << tol() << std::endl;
		if(isConverged) {
			np->unsetQuantities();
			if(y_inout != y.get()) assign( y_inout, *y );  // Assign the solution if we have to
			if(out) {
				*out << "\nWe have converged :-)\n"
					 << "\n||y||inf = " << norm_inf(*y) << std::endl;
				if(dumpAll) *out << "\ny =\n" << *y;
				*out << "\nExiting SimpleNewtonSolver::solve(...)\n";
			}
			return Solvers::SolveReturn(Solvers::SOLVED_TO_TOL,newtonIter);
		}
		if(out) *out << "\nWe have to keep going :-(\n";
		// Compute the Jacobian if we have not already
		if(newtonIter > 1) {
			if(out) *out << "\nComputing the Jacobian DcDy at current point ...\n";
			np->calc_DcDy(*y,NULL,false);
			if(out && dumpAll) {
				*out << "\nDcDy =\n" << *DcDy;
			}
		}
		// Compute the newton step: dy = -inv(DcDy)*c
		if(out) *out << "\nComputing the Newton step: dy = - inv(DcDy)*c ...\n";
		assign( dy.get(), 0.0 );         // Initial guess for the linear solve
		DcDy->solve(opDcDy,*c,dy.get()); // Solve: DcDy*dy = c
		Vt_S( dy.get(), -1.0 );          // dy *= -1.0
		if(out) *out << "\n||dy||inf = " << norm_inf(*dy) << std::endl;
		if(out && dumpAll) *out << "\ndy =\n" << *dy;
		// Perform backtracking armijo line search
		if(out) *out << "\nStarting backtracking line search iterations ...\n";
		const Scalar Dphi = -2.0*phi; // D(phi(y+alpha*dy))/D(alpha) at alpha=0.0 => -dot(c,c)/n: where dy = -inv(DcDy)*c
		Scalar alpha = 1.0; // Try a full step initially since it will eventually be accepted near solution
		int lineSearchIter;
		for( lineSearchIter = 1; lineSearchIter <= maxLineSearchIter(); ++lineSearchIter ) {
			if(out) *out << "\n*** lineSearchIter = " << lineSearchIter << std::endl;
			// y_new = y + alpha*dy
			assign( y_new.get(), *y ); Vp_StV( y_new.get(), alpha, *dy );
			if(out) *out << "\n||y_new||inf = " << norm_inf(*y_new) << std::endl;
			if(out && dumpAll) *out << "\ny_new =\n" << *y_new;
			// Compute the residual at the updated point
			np->calc_c(*y_new,NULL);
			if(out && dumpAll) *out << "\nc_new =\n" << *c;
			const Scalar phi_new = 0.5*dot(*c,*c)/n, phi_frac = phi + alpha * eta() * Dphi;
			const bool acceptPoint = (phi_new <= phi_frac);
			if(out) *out << "\nphi_new = 0.5*dot(c_new,c_new)/n = " << phi_new << std::endl
						 << "\nphi_new = " << phi_new << ( isConverged ? " <= " : " > " )
						 << "phi + alpha * eta * Dphi = " << phi << " + " << alpha << " * " << eta() << " * " << Dphi
						 << " = " << phi_frac << std::endl;
			if( acceptPoint ) {
				if(out) *out << "\nAccepting the current step with step length alpha = " << alpha << "!\n";
				break;
			}
			// Backtrack
			alpha *= 0.5;
		}
		// Check for line search failure
		if( lineSearchIter > maxLineSearchIter() ) {
			if(out) *out << "\nlineSearchIter = " << lineSearchIter << " > maxLineSearchIter = " << maxLineSearchIter()
						 << ": Terminating algorithm (throw SolverBreakDown)...\n";
			THROW_EXCEPTION(
				lineSearchIter > maxLineSearchIter(), Solvers::Exceptions::SolverBreakdown
				,"lineSearchIter = " << lineSearchIter << " > maxLineSearchIter = " << maxLineSearchIter()
				<< ": Terminating algorithm!" );
		}
		// Take the Newton step
		std::swap<mmp::ref_count_ptr<Vector<Scalar> > >( y_new, y ); // Swap y_new and y
	}
	np->unsetQuantities();
	// Failure!
	if(out) *out << "\nnewtonIter = " << newtonIter << " > maxNewtonIter = " << maxNewtonIter()
				 << ": Terminating algorithm ...\n";
	return Solvers::SolveReturn(Solvers::MAX_ITER_EXCEEDED,newtonIter);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP
