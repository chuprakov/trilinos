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
// TSFCoreSolversGMRESSolver.hpp

#ifndef TSFCORE_SOLVERS_GMRES_SOLVER_HPP
#define TSFCORE_SOLVERS_GMRES_SOLVER_HPP

//#define TSFCORE_GMRES_HACKED_PRINT_STATEMENTS

#include "TSFCoreSolversGMRESSolverDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreSolversConvergenceTester.hpp"
#include "TSFCoreTestingTools.hpp"
#include "Teuchos_BLAS.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
GMRESSolver<Scalar>::GMRESSolver(
	const Teuchos::RefCountPtr<std::ostream>   &out
	,bool                                      dump_all
	,int                                       default_max_iter
	,ScalarMagnitude                           breakdown_tol
	)
	:out_(out)
	,dump_all_(dump_all)
	,default_max_iter_(default_max_iter)
	,breakdown_tol_(breakdown_tol)
{}

// Overridden from SolverState

template<class Scalar>
Index GMRESSolver<Scalar>::totalNumSystems() const
{
	return totalNumSystems_;
}

template<class Scalar>
Index GMRESSolver<Scalar>::currNumSystems() const
{
	return 1; // Current implementation always solves one system at a time!
}

template<class Scalar>
int GMRESSolver<Scalar>::currIteration() const
{
	return curr_iter_;
}

template<class Scalar>
void GMRESSolver<Scalar>::currActiveSystems( Index activeSystems[] ) const
{
	activeSystems[0] = currActiveSystem_;
}

template<class Scalar>
void GMRESSolver<Scalar>::currEstRelResidualNorms( Scalar norms[] ) const
{
	norms[0] = rel_r_nrm_;
}

// Overridden from IterativeLinearSolver

template<class Scalar>
bool GMRESSolver<Scalar>::adjointRequired() const
{
	return false;
}

template<class Scalar>
SolveReturn GMRESSolver<Scalar>::solve(
	const LinearOp<Scalar> &M, const ETransp M_trans, const MultiVector<Scalar> &Y, MultiVector<Scalar> *X
	,const Scalar a, const int max_iter_in, ConvergenceTester<Scalar> *convTester
	,const LinearOp<Scalar> *M_tilde_left_inv, const ETransp M_tilde_left_inv_trans
	,const LinearOp<Scalar> *M_tilde_right_inv, const ETransp M_tilde_right_inv_trans
	) const
{
	TEST_FOR_EXCEPTION(
		M_tilde_right_inv!=NULL || M_tilde_right_inv != NULL, std::logic_error
		,"Error, we can not handle preconditioners yet!"
		);
	namespace mmp = MemMngPack;
	const VectorSpace<Scalar> &opM_domain     = ( M_trans == NOTRANS                ? *M.domain() : *M.range()  );
	const VectorSpace<Scalar> &opM_range      = ( M_trans == NOTRANS                ? *M.range()  : *M.domain() );
	//
	totalNumSystems_ = Y.domain()->dim();
	//
	// Validate input
	//
#ifdef _DEBUG
	const char func_name[] = "GMRESSolver<Scalar>::solve(...)";
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
		*get_out() << "\n*** Entering GMRESSolver<Scalar>::solve(...)\n" << std::setprecision(16);
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
	const int max_iter = ( max_iter_in == DEFAULT_MAX_ITER ? default_max_iter() : max_iter_in );
	//
	// Solve each linear system one at a time
	//
	bool all_solved = true;
	int max_iter_taken = 0;
	Teuchos::RefCountPtr<Vector<Scalar> >
		y = Y.range()->createMember(),
		x = Teuchos::null;
	r_ = opM_range.createMember();              // Used as workspace in localSolve
	V_ = opM_domain.createMembers(max_iter+1);  // Used as workspace in localSolve
	for( int k = 1; k <= totalNumSystems_; ++k ) {
		currActiveSystem_ = k;
		assign( &*y, *Y.col(k) );  Vt_S( &*y, a );  // y = a*Y.col(k)  (copy)
		x = X->col(k);                              // x = X.col(k)    (view)
		if(convTester) convTester->reset();
		const SolveReturn single_solve_return =
			this->localSolve(
				M                       // M
				,M_trans                // M_trans
				,*y                     // b
				,&*x                    // x
				,max_iter               // max_iter
				,convTester             // convTester
				);
		if( single_solve_return.solve_status == MAX_ITER_EXCEEDED ) all_solved = false;
		if( single_solve_return.num_iter > max_iter_taken ) max_iter_taken = single_solve_return.num_iter;
		const Scalar currEstRelResidualNorm = rel_r_nrm_;
		if(get_out().get()) {
			*get_out()
				<< "\nGMRESSolver<Scalar>::solve(...) : k = "<< k << " system : "
				<< ( single_solve_return.solve_status == SOLVED_TO_TOL ? "Solved for x(" : "Did not solve for x(" ) << k << ")"
				<< "\n  Number of GMRES iteratins taken          = " << single_solve_return.num_iter
				<< "\n  Final GMRES preconditioned residual norm = " << currEstRelResidualNorm << std::endl;
			if(dump_all()) {
				*get_out() << "\nX =\n" << *X;
			}
			*get_out() << "\n*** Leaving GMRESSolver<Scalar>::solve(...)\n";
		}
	}
	//
	// Return the solution
	//
	if(get_out().get()) {
		*get_out()
			<< "\nGMRESSolver<Scalar>::solve(...) : " << ( all_solved ? "Solved for X" : "Did not solve for X" );
		if(dump_all()) {
			*get_out() << "\nX =\n" << *X;
		}
		*get_out() << "\n*** Leaving GMRESSolver<Scalar>::solve(...)\n";
	}
	return SolveReturn( all_solved ? SOLVED_TO_TOL : MAX_ITER_EXCEEDED , max_iter_taken );
}

template<class Scalar>
Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> >
GMRESSolver<Scalar>::clone() const
{
	return Teuchos::rcp( new GMRESSolver<Scalar>(*this) );
}

// private

template<class Scalar>
SolveReturn GMRESSolver<Scalar>::localSolve(
	const LinearOp<Scalar>               &M
	,const ETransp                       M_trans
	,const Vector<Scalar>                &b
	,Vector<Scalar>                      *x
	,const int                           max_iter
	,ConvergenceTester<Scalar>           *convTester
	) const
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	const Scalar one = ST::one();
	// 
	// Check compatability of linear operator with rhs and solution vector
	//
	const VectorSpace<Scalar> &M_domain = ( ( M_trans == NOTRANS ) ? *M.domain() : *M.range()  );
	const VectorSpace<Scalar> &M_range  = ( ( M_trans == NOTRANS ) ? *M.range()  : *M.domain() );
	const bool
		domain_compatable = M_domain.isCompatible( *x->space() ),
		range_compatable  = M_range.isCompatible( *b.space() );
	TEST_FOR_EXCEPT(!(domain_compatable && range_compatable));
	//
	//  Initialize internal data structures
	//		
	curr_iter_ = 0;
	//
	if( max_iter+1 != H_.numRows() ) {
		H_.shapeUninitialized(max_iter+1,max_iter);
		z_.resize(max_iter+1);
		cs_.resize(max_iter);
		sn_.resize(max_iter);
	}
	//
	// Check if the RHS is zero
	//
#ifdef TSFCORE_GMRES_USE_DOT_FOR_SCALAR_PROD	
	const Scalar norm_b = norm_2( b );
#else
	const Scalar norm_b = norm( b );
#endif
#ifdef TSFCORE_GMRES_HACKED_PRINT_STATEMENTS
		std::cout << "\nGMRESSolver::solve(): ||b|| = " << norm_b << std::endl;
#endif
	bool isConverged = false;
	if( norm_b == ST::zero() ) {
		isConverged = true;
		assign( x, 0.0 );
	}
	if(!isConverged) {
		//
		// Determine the residual from the current solution: r = M*x - b 
		//
		assign( &*r_, b );
		const ScalarMagnitude norm_x = norm_inf(*x);
		if( norm_x == ST::zero() )
			Vt_S( &*r_, -one );
		else
			M.apply( M_trans, *x, &*r_, one, -one );
		//
		// Set up initial vector.
		//
#ifdef TSFCORE_GMRES_USE_DOT_FOR_SCALAR_PROD	
		r0_nrm_ = norm_2( *r_ );
#else
		r0_nrm_ = norm( *r_ );
#endif
#ifdef TSFCORE_GMRES_HACKED_PRINT_STATEMENTS
		std::cout << "\nGMRESSolver::solve(): ||r0|| = " << r0_nrm_ << std::endl;
#endif
		rel_r_nrm_ = r0_nrm_ / norm_b;
		//
		// Do the GMRES iterations
		//
		bool solverBreakdown = false;
		while( !isConverged ) {
			//
			// Convergence check
			//
#ifdef TSFCORE_GMRES_HACKED_PRINT_STATEMENTS
			std::cout << "\nGMRESSolver::solve(): ||r("<<curr_iter_<<")||/(1+||b||) = " << rel_r_nrm_ << std::endl;
#endif
			if(convTester) {
				bool isConvergedArg[1];
				convTester->convStatus(*this,1,isConvergedArg);
				if(isConvergedArg[0]) {
					isConverged = true;
					break;
				}
			}
			if(solverBreakdown) {
				isConverged = true;
				break;
			}
			if(curr_iter_ == max_iter)
				break;
			//
			// Setup for first GMRES iteration
			//
			if(curr_iter_==0) {
				z_[0] = r0_nrm_;
				Teuchos::RefCountPtr<Vector<Scalar> >
					w = V_->col(1);		             // Get a mutable view of the first column of V_.
				assign( &*w, *r_ );              // Copy r to the first column of V_.
				Vt_S( &*w, one/r0_nrm_ );        // v_1 = r_0 / ||r_0||
				w = Teuchos::null;               // Forces V_ to be modified
			}
			//
			// Do the iteration
			//
			solverBreakdown = doIteration( M, M_trans );
			//
			// Increment the iteration counter to the iteration number just completed
			//
			curr_iter_++;
		}
		//
		// Solve least squares problem.
		//
		Teuchos::BLAS<int, Scalar> blas;
		blas.TRSM(
			Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, 
			curr_iter_, 1, one, H_.values(), max_iter+1, &z_[0], max_iter+1
			);
		//
		// Compute the new solution.
		//
		if( curr_iter_ > 0 ) {
			Teuchos::RefCountPtr<MultiVector<Scalar> > V = V_->subView(Range1D(1,curr_iter_));
			SerialVector<Scalar> z_vec( Teuchos::rcp(&z_[0],false), 1, curr_iter_ );
			V->apply( NOTRANS, z_vec, x, -one, one );
		}
	}
	return SolveReturn( curr_iter_ >= max_iter ? MAX_ITER_EXCEEDED : SOLVED_TO_TOL, curr_iter_ );
}

template<class Scalar>
bool GMRESSolver<Scalar>::doIteration( const LinearOp<Scalar> &M, const ETransp M_trans ) const
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	const Scalar one = ST::one();
	int i;
	Teuchos::BLAS<int, Scalar> blas;
	Teuchos::SerialDenseMatrix<int, Scalar> &H = H_;
	std::valarray<Scalar> &z = z_, &cs = cs_, &sn = sn_;
	// 
	Teuchos::RefCountPtr<Vector<Scalar> >
		w = V_->col(curr_iter_+2);                                           // w = v_{j+1}
	M.apply( M_trans, *V_->col(curr_iter_+1), &*w );                       // w = M * v_{j}	
	//
	// Perform MGS to orthogonalize new Krylov vector.
	//
	for( i = 0; i < curr_iter_ + 1; i++ ) {
#ifdef TSFCORE_GMRES_USE_DOT_FOR_SCALAR_PROD	
		H( i, curr_iter_ ) = dot( *w, *V_->col(i+1) );                       // h_{i,j} = ( w, v_{i} )
#else
		H( i, curr_iter_ ) = w->space()->scalarProd( *w, *V_->col(i+1) );    // h_{i,j} = ( w, v_{i} )
#endif
		Vp_StV( &*w, -H( i, curr_iter_ ), *V_->col(i+1) );                   // w = w - h_{i,j} * v_{i}
	}
	const Scalar H_jp1_j = H( curr_iter_+1, curr_iter_ ) =
#ifdef TSFCORE_GMRES_USE_DOT_FOR_SCALAR_PROD	
		norm_2( *w );                                                        // h_{j+1,j} = || w ||
#else
		norm( *w );                                                          // h_{j+1,j} = || w ||
#endif
	const bool breakdown = ST::magnitude(H_jp1_j) <= breakdown_tol();
	if(!breakdown) Vt_S( &*w, one / H_jp1_j );                             // v_{j+1} = w / h_{j+1,j}			
	//
	// Apply previous Givens rotations
	//
	for( i = 0; i < curr_iter_; i++ ) {
		const Scalar temp = cs[i]*H( i, curr_iter_ ) + sn[i]*H( i+1, curr_iter_ );
		H( i+1, curr_iter_ ) = -sn[i]*H( i, curr_iter_ ) + cs[i]*H( i+1, curr_iter_ );
		H( i, curr_iter_ ) = temp;
	}
	//
	// Calculate new Givens rotation
	//
	blas.ROTG(
		&H( curr_iter_, curr_iter_ ), &H( curr_iter_+1, curr_iter_ ), 
		&cs[curr_iter_], &sn[curr_iter_]
		);
	//
	// Update RHS and residual w/ new transform and compute residual norm.
	//
	z[curr_iter_+1] = -sn[curr_iter_]*z[curr_iter_];
	z[curr_iter_] *= cs[curr_iter_];
#ifdef TSFCORE_GMRES_HACKED_PRINT_STATEMENTS
	std::cout << "\nGMRESSolver::solve(): z["<<curr_iter_+1<<"] = " << z[curr_iter_+1] << std::endl;
#endif
	rel_r_nrm_ = ST::magnitude( z[curr_iter_+1] ) / r0_nrm_;
	//
	// Return that the solver has not broken down
	//
	return breakdown;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_GMRES_SOLVER_HPP
