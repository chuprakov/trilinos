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

// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreSolversSimpleGMRESSolver.hpp

#ifndef TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP
#define TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP

#include "TSFCoreSolversSimpleGMRESSolverDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreSolversConvergenceTester.hpp"
#include "TSFCoreSolversNormedConvergenceTester.hpp"
#include "TSFCoreTestingTools.hpp"
#include "check_nan_inf.h"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
SimpleGMRESSolver<Scalar>::SimpleGMRESSolver(
	const Teuchos::RefCountPtr<std::ostream>   &out
	,bool                                      dump_all
	,int                                       default_max_iter
	,Scalar                                    default_tol
	)
	:out_(out)
	,dump_all_(dump_all)
	,default_max_iter_(default_max_iter)
	,default_tol_(default_tol)
{}

// Overridden from SolverState

template<class Scalar>
Index SimpleGMRESSolver<Scalar>::totalNumSystems() const
{
	assert(0);
	return 0;
}

template<class Scalar>
Index SimpleGMRESSolver<Scalar>::currNumSystems() const
{
	assert(0);
	return 0;
}

template<class Scalar>
int SimpleGMRESSolver<Scalar>::currIteration() const
{
	assert(0);
	return 0;
}

template<class Scalar>
void SimpleGMRESSolver<Scalar>::currActiveSystems( Index activeSystems[] ) const
{
	assert(0);
}

template<class Scalar>
void SimpleGMRESSolver<Scalar>::currEstRelResidualNorms( Scalar norms[] ) const
{
	assert(0);
}

// Overridden from IterativeLinearSolver

template<class Scalar>
bool SimpleGMRESSolver<Scalar>::adjointRequired() const
{
	return false;
}

template<class Scalar>
SolveReturn SimpleGMRESSolver<Scalar>::solve(
	const LinearOp<Scalar> &M, ETransp M_trans, const MultiVector<Scalar> &Y, MultiVector<Scalar> *X
	,const Scalar a, const int max_iter_in, ConvergenceTester<Scalar> *convTester
	,const LinearOp<Scalar> *M_tilde_left_inv, const ETransp M_tilde_left_inv_trans
	,const LinearOp<Scalar> *M_tilde_right_inv, ETransp M_tilde_right_inv_trans
	) const
{
	TEST_FOR_EXCEPTION(
		M_tilde_right_inv!=NULL || M_tilde_right_inv != NULL, std::logic_error
		,"Error, we can not handle preconditioners yet!"
		);
	namespace mmp = MemMngPack;
	const VectorSpace<Scalar> &opM_domain     = ( M_trans == NOTRANS                ? *M.domain() : *M.range()  );
	const VectorSpace<Scalar> &opM_range      = ( M_trans == NOTRANS                ? *M.range()  : *M.domain() );
	const ETransp      opM_notrans            = ( M_trans == NOTRANS                ? NOTRANS     : TRANS       );
	//const ETransp      opM_trans              = ( M_trans == NOTRANS                ? TRANS       : NOTRANS     );
	//const ETransp      opM_tilde_inv_notrans  = ( M_tilde_left_inv_trans == NOTRANS ? NOTRANS     : TRANS       );
	//const ETransp      opM_tilde_inv_trans    = ( M_tilde_left_inv_trans == NOTRANS ? TRANS       : NOTRANS     );
	//
	TEST_FOR_EXCEPTION(
		opM_notrans!=NOTRANS, std::logic_error
		,"Error, we can not handle the transpose argument yet!";
		);
	//
	const int totalNumSystems = Y.domain()->dim();
	//
	// Validate input
	//
#ifdef _DEBUG
	const char func_name[] = "SimpleGMRESSolver<Scalar>::solve(...)";
	TEST_FOR_EXCEPTION(X==NULL,std::invalid_argument,": Error!");
	if(M_tilde_left_inv) {
		const VectorSpace<Scalar> &opM_tilde_inv_domain = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->domain() : *M_tilde_left_inv->range()  );
		const VectorSpace<Scalar> &opM_tilde_inv_range  = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->range()  : *M_tilde_left_inv->domain() );
		const bool
			domain_compatible = opM_domain.isCompatible(opM_tilde_inv_range),
			range_compatible =  opM_range.isCompatible(opM_tilde_inv_domain);
		TEST_FOR_EXCEPTION(
			!(domain_compatible && range_compatible), TSFCore::Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_inv) do are not compatible!");
	}
	bool is_compatible = opM_domain.isCompatible(*X->range());
	TEST_FOR_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the op(M).domain() not compatible with X->range()!");
	is_compatible = opM_range.isCompatible(*Y.range());
	TEST_FOR_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the op(M).range() not compatible with Y.range()!");
	is_compatible = X->domain()->isCompatible(*Y.domain());
	TEST_FOR_EXCEPTION(
		!is_compatible, TSFCore::Exceptions::IncompatibleVectorSpaces
		,func_name<<": Error, the X->domain() not compatible with Y.domain()!");
#endif
	if(get_out().get()) {
		*get_out() << "\n*** Entering SimpleGMRESSolver<Scalar>::solve(...)\n" << std::setprecision(16);
		if(dump_all()) {
			*get_out() << "\nM =\n" << M;
			*get_out() << "\nM_trans = " << toString(M_trans) << std::endl;
			*get_out() << "\nY =\n" << Y;
			*get_out() << "\nX =\n" << *X;
			*get_out() << "\na = " << a << std::endl;
		}
	}
	//
	// Resolve default parameters
	//
	NormedConvergenceTester<Scalar>
		*normedConvTester = dynamic_cast<NormedConvergenceTester<Scalar>*>(convTester);
	const int max_iter = ( max_iter_in == DEFAULT_MAX_ITER ? default_max_iter() : max_iter_in );
	const Scalar tol = ( normedConvTester ? normedConvTester->tol() : default_tol_ );
	//
	// Solve each linear system one at a time
	//
	bool all_solved = true;
	int max_iter_taken = 0;
	Scalar minMaxErr = Scalar(1e+50);
	Teuchos::RefCountPtr<Vector<Scalar> >
		y = Y.range()->createMember(),
		x = Teuchos::null;
	for( int k = 1; k <= totalNumSystems; ++k ) {
		assign( &*y, *Y.col(k) );  Vt_S( &*y, a );  // y = a*Y.col(k)  (copy)
		x = X->col(k);                              // x = X.col(k)    (view)
		const SolveReturn single_solve_return
			= solver_.solve(
				M                       // Op
				,*y                     // b
				,&*x                    // curr_soln
				,M_trans                // Op_trans
				,max_iter               // max_iter_in
				,tol                    // tol_in
				);
		if( single_solve_return.solve_status == MAX_ITER_EXCEEDED ) all_solved = false;
		if( single_solve_return.num_iter > max_iter_taken ) max_iter_taken = single_solve_return.num_iter;
		const Scalar currEstRelResidualNorm = solver_.currEstRelResidualNorm();
		if(currEstRelResidualNorm < minMaxErr) minMaxErr = currEstRelResidualNorm;
		if(get_out().get()) {
			*get_out()
				<< "\nSimpleGMRESSolver<Scalar>::solve(...) : k = "<< k << " system : "
				<< ( single_solve_return.solve_status == SOLVED_TO_TOL ? "Solved for x(" : "Did not solve for x(" ) << k << ")"
				<< "\n  Number of GMRES iteratins taken          = " << single_solve_return.num_iter
				<< "\n  Final GMRES preconditioned residual norm = " << currEstRelResidualNorm << std::endl;
			if(dump_all()) {
				*get_out() << "\nX =\n" << *X;
			}
			*get_out() << "\n*** Leaving SimpleGMRESSolver<Scalar>::solve(...)\n";
		}
	}
	if(normedConvTester) normedConvTester->minMaxErr(minMaxErr);
	//
	// Return the solution
	//
	if(get_out().get()) {
		*get_out()
			<< "\nSimpleGMRESSolver<Scalar>::solve(...) : " << ( all_solved ? "Solved for X" : "Did not solve for X" );
		if(dump_all()) {
			*get_out() << "\nX =\n" << *X;
		}
		*get_out() << "\n*** Leaving SimpleGMRESSolver<Scalar>::solve(...)\n";
	}
	return SolveReturn( all_solved ? SOLVED_TO_TOL : MAX_ITER_EXCEEDED , max_iter_taken );
}

template<class Scalar>
Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> >
SimpleGMRESSolver<Scalar>::clone() const
{
	return Teuchos::rcp( new SimpleGMRESSolver<Scalar>(*this) );
}


} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_SIMPLE_GMRES_SOLVER_HPP
