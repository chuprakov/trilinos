// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// //////////////////////////////////////////////////////////////////
// TSFCoreNonlinSimpleNewtonSolver.hpp

#ifndef TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP
#define TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP

#include "TSFCoreNonlinSimpleNewtonSolverDecl.hpp"
#include "TSFCoreTestingTools.hpp"
#include "TSFCoreSolversSummaryOutputter.hpp"
#include "TSFCoreNonlinNonlinearProblem.hpp"
#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "AbstractFactory.hpp"

namespace TSFCore {
namespace Nonlin {

template<class T>
void my_sad_swap( T &t1, T &t2 ) {
	// RAB: 05/03/2004: Replaced HKT's incode swap to call this function.
	// ToDo: Once std::swap<> is supported get rid of this function.
	T t_tmp;
	t_tmp = t1;
	t1 = t2;
	t2 = t_tmp;
}


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
SimpleNewtonSolver<Scalar>::solve(
	NonlinearProblemFirstOrder<Scalar> *np
	,Vector<Scalar> *y_inout, std::ostream *out, bool dumpAll
	) const
{
	typedef Teuchos::ScalarTraits<Scalar>  ST;
	const std::string &leadstr = leadingSummaryOutputStr();
	// Setup summary outputting
	Teuchos::RefCountPtr<Solvers::SummaryOutputter<Scalar> > souter;
	std::ostream *sout = get_summaryOut().get();
	if(sout) {
		souter = Teuchos::rcp(
			new Solvers::SummaryOutputter<Scalar>(
				get_summaryOut()
				,leadstr+std::string("  ")
				)
			);
	}
	if(out) *out << "\n*** Entering SimpleNewtonSolver::solve(...) ...\n";
 	if(sout) *sout << leadstr << "newton_iter=0: Starting nonlinear solve ...\n";
	// Create the data for the problem
	const VectorSpace<Scalar>                         &space_c = *np->space_c();
	const VectorSpace<Scalar>                         &space_y = *np->space_y();
	Teuchos::RefCountPtr<Vector<Scalar> >             y        = Teuchos::rcp(y_inout,false);   // Initial guess/solution
	Teuchos::RefCountPtr<Vector<Scalar> >             c        = space_c.createMember();        // Residual
	Teuchos::RefCountPtr<LinearOpWithSolve<Scalar> >  DcDy     = np->factory_DcDy()->create();  // Jacobian object
	ETransp                                           opDcDy   = np->opDcDy();                  // Transpose argument for DcDy
	Teuchos::RefCountPtr<Vector<Scalar> >             dy       = space_y.createMember();        // Newton step for y
	Teuchos::RefCountPtr<Vector<Scalar> >             y_new    = space_y.createMember();        // Trial point for y
	// Compute the initial starting point
	np->unsetQuantities(); np->set_c(&*c); np->set_DcDy(&*DcDy); // These pointers will be maintained throughout
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
	try {
		for( newtonIter = 1; newtonIter <= maxNewtonIter(); ++newtonIter ) {
			if(out) *out << "\n*** newtonIter = " << newtonIter << std::endl;
			// Check convergence
			if(out) *out << "\nChecking for convergence ... : ";
			const Scalar phi = space_c.scalarProd(*c,*c), sqrt_phi = sqrt(phi); // merit function: phi(c) = <c,c>
			const bool isConverged = sqrt_phi <= tol();
			if(out) *out << "sqrt(phi) = sqrt(<c,c>) = ||c|| = " << sqrt_phi << ( isConverged ? " <= " : " > " ) << "tol = " << tol() << std::endl;
			if(sout) *sout
				<< leadstr << "newton_iter="<<newtonIter<<": Check convergence: ||c|| = "
				<< sqrt_phi << ( isConverged ? " <= " : " > " ) << "tol = " << tol() << ( isConverged ? ", Converged!!!" : "" ) << std::endl;
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
			if(sout) *sout	<< leadstr << "newton_iter="<<newtonIter<<": Computing Newton step ...\n";
			assign( &*dy, 0.0 );                      // Initial guess for the linear solve
			DcDy->solve(opDcDy,*c,&*dy,souter.get()); // Solve: DcDy*dy = c
			Vt_S( &*dy, -1.0 );                       // dy *= -1.0
			if(out) *out << "\n||dy||inf = " << norm_inf(*dy) << std::endl;
			if(out && dumpAll) *out << "\ndy =\n" << *dy;
			// Perform backtracking armijo line search
			if(out) *out << "\nStarting backtracking line search iterations ...\n";
			if(sout) *sout	<< leadstr << "newton_iter="<<newtonIter<<": Starting backtracking line search ...\n";
			const Scalar Dphi = -2.0*phi; // D(phi(y+alpha*dy))/D(alpha) at alpha=0.0 => -2.0*<c,c>: where dy = -inv(DcDy)*c
			Scalar alpha = 1.0; // Try a full step initially since it will eventually be accepted near solution
			int lineSearchIter;
			for( lineSearchIter = 1; lineSearchIter <= maxLineSearchIter(); ++lineSearchIter ) {
				if(out) *out << "\n*** lineSearchIter = " << lineSearchIter << std::endl;
				// y_new = y + alpha*dy
				assign( &*y_new, *y ); Vp_StV( &*y_new, alpha, *dy );
				if(out) *out << "\n||y_new||inf = " << norm_inf(*y_new) << std::endl;
				if(out && dumpAll) *out << "\ny_new =\n" << *y_new;
				// Compute the residual at the updated point
				np->calc_c(*y_new,NULL);
				if(out && dumpAll) *out << "\nc_new =\n" << *c;
				const Scalar phi_new = space_c.scalarProd(*c,*c), phi_frac = phi + alpha * eta() * Dphi;
				if(out) *out << "\nphi_new = <c_new,c_new> = " << phi_new << std::endl;
				if( Teuchos::ScalarTraits<Scalar>::isnaninf(phi_new) ) {
					if(out) *out << "\nphi_new is not a valid number, backtracking (alpha = 0.1*alpha) ...\n";
					alpha *= 0.1;
					continue;
				}
				const bool acceptPoint = (phi_new <= phi_frac);
				if(out) *out
					<< "\nphi_new = " << phi_new << ( acceptPoint ? " <= " : " > " )
					<< "phi + alpha * eta * Dphi = " << phi << " + " << alpha << " * " << eta() << " * " << Dphi
					<< " = " << phi_frac << std::endl;
				if(sout) *sout
					<< leadstr << "newton_iter="<<newtonIter<<", ls_iter="<<lineSearchIter<<" : "
					<< "phi(alpha="<<alpha<<") = "<<phi_new<<(acceptPoint?" <=":" >")<<" armijo_cord = " << phi_frac << std::endl;
				if( acceptPoint ) {
					if(out) *out << "\nAccepting the current step with step length alpha = " << alpha << "!\n";
					break;
				}
				if(out) *out << "\nBacktracking (alpha = 0.5*alpha) ...\n";
				alpha *= 0.5;
			}
			// Check for line search failure
			if( lineSearchIter > maxLineSearchIter() ) {
				if(out) *out
					<< "\nlineSearchIter = " << lineSearchIter << " > maxLineSearchIter = " << maxLineSearchIter()
					<< ": Terminating algorithm (throw SolverBreakDown)...\n";
				if(sout) *sout
					<< std::endl << leadstr << "ls_iter = " << lineSearchIter << " > maxLineSearchIter = " << maxLineSearchIter()
					<< ": Terminating algorithm (throw SolverBreakDown)...\n";
				if(y_inout != y.get()) assign( y_inout, *y ); // Assign the final point
				TEST_FOR_EXCEPTION(
					lineSearchIter > maxLineSearchIter(), Solvers::Exceptions::SolverBreakdown
					,"lineSearchIter = " << lineSearchIter << " > maxLineSearchIter = " << maxLineSearchIter()
					<< ": Terminating algorithm!" );
			}
			// Take the Newton step
			my_sad_swap<Teuchos::RefCountPtr<Vector<Scalar> > >( y_new, y ); // Now y is the accepted point!
		}
	}
	catch (...) {
		if(out)
			*out << "\nSimpleNewtonSolver<Scalar>::solve(...): An exception was thrown, setting the last accepted point before rethrowing ...\n";
		if(y_inout != y.get()) assign( y_inout, *y ); // Assign the final point
		throw;
	}
	np->unsetQuantities();
	// Failure!
	if(out) *out << "\nnewtonIter = " << newtonIter << " > maxNewtonIter = " << maxNewtonIter()
				 << ": Terminating algorithm ...\n";
	if(y_inout != y.get()) assign( y_inout, *y ); // Assign the final point
	return Solvers::SolveReturn(Solvers::MAX_ITER_EXCEEDED,newtonIter);
}

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_SIMPLE_NEWTON_SOLVER_HPP
