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
// TSFCoreSolversBiCGSolver.hpp

#ifndef TSFCORE_SOLVERS_BICG_SOLVER_HPP
#define TSFCORE_SOLVERS_BICG_SOLVER_HPP

#include "TSFCoreSolversBiCGSolverDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "TSFCoreSolversConvergenceTester.hpp"
#include "TSFCoreTestingTools.hpp"
#include "check_nan_inf.h"
#include "Teuchos_TestForException.hpp"

namespace {

template<class Scalar>
void compressMember( TSFCore::Index k, TSFCore::Index k_curr, const Teuchos::RefCountPtr<TSFCore::MultiVector<Scalar> >& V )
{
	TSFCore::assign<Scalar>( V->col(k_curr+1).get(), *V->col(k+1) ); 
}

template<class T>
void compressMember( TSFCore::Index k, TSFCore::Index k_curr, std::vector<T>* a )
{
	(*a)[k_curr] = (*a)[k];
}

} // namespace

namespace TSFCore {
namespace Solvers {

template<class Scalar>
BiCGSolver<Scalar>::BiCGSolver(
	const out_ptr_t   &out
	,bool             dump_all
	,int              default_max_iter
	)
	:out_(out)
	,dump_all_(dump_all)
	,default_max_iter_(default_max_iter)
{}

// Overridden from SolverState

template<class Scalar>
Index BiCGSolver<Scalar>::totalNumSystems() const
{
	return totalNumSystems_;
}

template<class Scalar>
Index BiCGSolver<Scalar>::currNumSystems() const
{
	return currNumSystems_;
}

template<class Scalar>
int BiCGSolver<Scalar>::currIteration() const
{
	return currIteration_ - 1; // We check convergence before the iteration
}

template<class Scalar>
void BiCGSolver<Scalar>::currActiveSystems( Index activeSystems[] ) const
{
	std::copy( activeSystems_.begin(), activeSystems_.end(), activeSystems );
}

template<class Scalar>
void BiCGSolver<Scalar>::currEstRelResidualNorms( Scalar norms[] ) const
{
	if(!norms_updated_) {
		TSFCore::norms( *R_, &norms_[0] );
		for(int j=0;j<currNumSystems_;++j)
			norms_[j] /= rel_err_denom_[j];
		norms_updated_ = true;
		if(get_out().get()) {
			*get_out() << "\n||R|| =\n"; for(int j=0;j<currNumSystems_;++j) *get_out() << " " << norms_[j]; *get_out() << std::endl;
		}
	}
	std::copy( norms_.begin(), norms_.end(), norms );
}

// Overridden from IterativeLinearSolver

template<class Scalar>
bool BiCGSolver<Scalar>::adjointRequired() const
{
	return true;
}

//
// Here we will transform the problem from
//
//    1/a*op(M)*X[j] - Y[j] == 0
//
// to
//
//    Y_hat[j] - op(M_hat)*X_hat[j] == 0
//
// with preconditioner M_tilde_inv
//
//    where:
//        if M_tilde_right_inv == NULL
//            Y_hat       = a*Y
//            M_hat       = M
//            X_hat       = X
//            M_tilde_inv = M_tilde_left_inv
//        elif M_tilde_left_inv == NULL && M_tilde_right_inv != NULL
//            Y_hat       = a*Y
//            M_hat       = M*M_tilde_right_inv
//            X_hat       = inv(M_tilde_right_inv)*X
//            M_tilde_inv = NULL
//        elif M_tilde_left_inv != NULL && M_tilde_right_inv != NULL
//            Y_hat       = a*M_tilde_left_inv*Y
//            M_hat       = M_tilde_left_inv*M*M_tilde_right_inv
//            X_hat       = inv(M_tilde_right_inv)*X
//            M_tilde_inv = NULL
//        endif
//
// and then apply the BiCG method as exactly described
// in "Templates for the Solution of Linear Systems" except
// this is a simple multi-vector version.
//
template<class Scalar>
SolveReturn BiCGSolver<Scalar>::solve(
	const LinearOp<Scalar> &M, ETransp M_trans, const MultiVector<Scalar> &Y, MultiVector<Scalar> *X
	,const Scalar a, const int max_iter_in, ConvergenceTester<Scalar> *convTester
	,const LinearOp<Scalar> *M_tilde_left_inv, const ETransp M_tilde_left_inv_trans
	,const LinearOp<Scalar> *M_tilde_right_inv, ETransp M_tilde_right_inv_trans
	) const
{
	assert( M_tilde_right_inv == NULL );
	namespace mmp = MemMngPack;
	const VectorSpace<Scalar> &opM_domain     = ( M_trans == NOTRANS                ? *M.domain() : *M.range()  );
	const VectorSpace<Scalar> &opM_range      = ( M_trans == NOTRANS                ? *M.range()  : *M.domain() );
	const ETransp      opM_notrans            = ( M_trans == NOTRANS                ? NOTRANS     : TRANS       );
	const ETransp      opM_trans              = ( M_trans == NOTRANS                ? TRANS       : NOTRANS     );
	const ETransp      opM_tilde_inv_notrans  = ( M_tilde_left_inv_trans == NOTRANS ? NOTRANS     : TRANS       );
	const ETransp      opM_tilde_inv_trans    = ( M_tilde_left_inv_trans == NOTRANS ? TRANS       : NOTRANS     );
	totalNumSystems_ = Y.domain()->dim();
	//
	int j;
	//
	// Validate input
	//
#ifdef _DEBUG
	const char func_name[] = "BiCGSolver<Scalar>::solve(...)";
	TEST_FOR_EXCEPTION(X==NULL,std::invalid_argument,": Error!");
	bool adjoint_supported = M.opSupported(not_trans(M_trans));
	TEST_FOR_EXCEPTION(!adjoint_supported,std::invalid_argument,func_name<<": Non-transposed and transposed ops of M must be supported!");
	if(M_tilde_left_inv) {
		const VectorSpace<Scalar> &opM_tilde_inv_domain = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->domain() : *M_tilde_left_inv->range()  );
		const VectorSpace<Scalar> &opM_tilde_inv_range  = ( M_tilde_left_inv_trans==NOTRANS ? *M_tilde_left_inv->range()  : *M_tilde_left_inv->domain() );
		const bool
			domain_compatible = opM_domain.isCompatible(opM_tilde_inv_range),
			range_compatible =  opM_range.isCompatible(opM_tilde_inv_domain);
		TEST_FOR_EXCEPTION(
			!(domain_compatible && range_compatible), TSFCore::Exceptions::IncompatibleVectorSpaces
			,func_name<<": Error, the range and/or domain spaces of op(M) and op(M_tilde_inv) do are not compatible!");
		bool adjoint_supported = M_tilde_left_inv->opSupported(not_trans(M_tilde_left_inv_trans));
		TEST_FOR_EXCEPTION(!adjoint_supported,std::invalid_argument,func_name<<": Non-transposed and transposed ops of M_tilde_inv must be supported!");
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
		*get_out() << "\n*** Entering BiCGSolver<Scalar>::solve(...)\n" << std::setprecision(16);
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
	// Setup storage and initialize the algorithm
	//
	currNumSystems_ = totalNumSystems_;
	activeSystems_.resize(currNumSystems_); for(j=0;j<currNumSystems_;++j) activeSystems_[j] = j+1;
	X_curr_ = X->clone_mv();
	R_ = opM_range.createMembers(currNumSystems_);   R_tilde_ = opM_domain.createMembers(currNumSystems_);
	Q_ = opM_range.createMembers(currNumSystems_);   Q_tilde_ = opM_domain.createMembers(currNumSystems_);
	Z_ = opM_domain.createMembers(currNumSystems_);  Z_tilde_ = opM_range.createMembers(currNumSystems_);
	P_ = opM_domain.createMembers(currNumSystems_);  P_tilde_ = opM_range.createMembers(currNumSystems_);
	rho_.resize(currNumSystems_);   rho_old_.resize(currNumSystems_); beta_.resize(currNumSystems_);
	gamma_.resize(currNumSystems_); alpha_.resize(currNumSystems_);   norms_.resize(currNumSystems_);	
	// R^{0} = a*Y
	assign( R_.get(), 0.0 ); update( a, Y, R_.get() );
	// Denominator for relative error
	rel_err_denom_.resize(currNumSystems_);
	norms( *R_, &rel_err_denom_[0] );
	for(j=0;j<currNumSystems_;++j) rel_err_denom_[j] += 1.0;
	// R^{0} += - op(M)*X^{0}
	M.apply(M_trans,*X,R_.get(),-1.0,1.0);
	// Choose R_tilde^{0}
	randomize( -1.0, 1.0, R_tilde_.get() );  // R_tilde^{0} = rand(-1,+1) (what else can we do?)
	//
	// Perform the iterations
	//
	std::valarray<bool> isConverged(totalNumSystems_); // Note: std::vector<bool>::operator[]() does not return &bool
	for( currIteration_ = 1; currIteration_ <= max_iter; ++currIteration_ ) {
		// Check convergence
		if(convTester) {
			if(get_out().get()) *get_out() << "\nChecking convergence (iteration " << currIteration_ << ") ...";
			// Convergence test
			norms_updated_ = false;
			convTester->convStatus(*this,currNumSystems_,&isConverged[0]);
			// Compress
			compress(&isConverged[0],X);
			// Fully converged?
			bool not_converged = false;
			for(j=0;j<currNumSystems_;++j) { if(!isConverged[j]) not_converged = true; }
			if(!not_converged) {
				if(get_out().get()) {
					*get_out() << ": Converged on iteration " << currIteration_ << "!\n";
					if(dump_all()) *get_out() << "\nX =\n" << *X;
					*get_out() << "\nexiting BiCGSolver<Scalar>::solve(...)\n";
				}
				return SolveReturn(SOLVED_TO_TOL,currIteration_);  // We have converged!
			}
			else {
				if(get_out().get()) *get_out() << ": Not converged!\n";
			}
		}
		if( currIteration_ > 1 ) { for(j=0;j<currNumSystems_;++j) rho_old_[j] = rho_[j]; } // rho_{i-1}[j] = rho_{i-2}[j]
		// Do the Iteration
		Teuchos::RefCountPtr<MultiVector<Scalar> > X_curr = X_curr_->subView(Range1D(1,currNumSystems_));
		doIteration(M,opM_notrans,opM_trans,X_curr.get(),a,M_tilde_left_inv,opM_tilde_inv_notrans,opM_tilde_inv_trans);
	}
	return SolveReturn(MAX_ITER_EXCEEDED,currIteration_);
}

template<class Scalar>
Teuchos::RefCountPtr<const IterativeLinearSolver<Scalar> >
BiCGSolver<Scalar>::clone() const
{
	return Teuchos::rcp( new BiCGSolver<Scalar>(*this) );
}

// private

template<class Scalar>
void BiCGSolver<Scalar>::doIteration(
	const LinearOp<Scalar> &M, ETransp opM_notrans, ETransp opM_trans, MultiVector<Scalar> *X, Scalar a
	,const LinearOp<Scalar> *M_tilde_inv, ETransp opM_tilde_inv_notrans, ETransp opM_tilde_inv_trans
	) const
{
	if(get_out().get()) {
		*get_out() << "\n*** currIteration = " << currIteration_ << std::endl;
		if(dump_all()) {
			*get_out() << "\nX =\n" << *X;
			*get_out() << "\nR =\n" << *R_;
			*get_out() << "\nR_tilde =\n" << *R_tilde_;
		}
	}
	const Index m = currNumSystems_;
	int j;
	if( M_tilde_inv ) { // Preconditioner is available
		M_tilde_inv->apply( opM_tilde_inv_notrans, *R_,       Z_.get()       );  // M_tilde_inv*R^{i-1}              -> Z^{i-1}
		M_tilde_inv->apply( opM_tilde_inv_trans,   *R_tilde_, Z_tilde_.get() );  // M_tilde_inv'R_tilde^{i-1}        -> Z_tilde^{i-1}
	}
	else {          // No preconditioner is available
		assign( Z_.get(),       *R_        );                                    // R^{i-1}                          -> Z^{i-1}
		assign( Z_tilde_.get(), *R_tilde_  );                                    // R_tilde^{i-1}                    -> Z_tilde^{i-1}
	}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nZ =\n" << *Z_;
		*get_out() << "\nZ_tilde =\n" << *Z_tilde_;
	}
	dots( *Z_, *R_tilde_, &rho_[0] );                        // rho_{i-1}[j] = Z^{i-1}[j]' * R_tilde^{i-1}[j]
	if(get_out().get() && dump_all()) {
		*get_out() << "\nrho =\n"; for(j=0;j<m;++j) *get_out() << " " << rho_[j]; *get_out() << std::endl;
	}
	for(j=0;j<m;++j) { // Check for failure: rho_{i-1} = 0
		TEST_FOR_EXCEPTION(
			rho_[j] == 0.0, Exceptions::SolverBreakdown
			,"BiCGSolver<Scalar>::solve(...): Error, rho["<<j<<"] = 0.0, the method has failed!"
			);
	}
	if( currIteration_ == 1 ) {
		assign( P_.get(),       *Z_       );                  // Z^{i-1}                                         -> P^{i}
		assign( P_tilde_.get(), *Z_tilde_ );                  // Z_tilde^{i-1}                                   -> P_tilde^{i}
	}
	else {
		for(j=0;j<m;++j) beta_[j] = rho_[j]/rho_old_[j];      // rho_{i-1}[j]/rho_{i-2}[j]                        -> beta_{i-1}[j]
		update( *Z_,       &beta_[0], 1.0, P_.get()       );  // Z^{i-1}[j]       + beta_{i-1} * P^{i-1}[j]       -> P^{i}[j]
		update( *Z_tilde_, &beta_[0], 1.0, P_tilde_.get() );  // Z_tilde^{i-1}[j] + beta_{i-1} * P_tilde^{i-1}[j] -> P_tilde^{i}[j]
	}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nP =\n" << *P_;
		*get_out() << "\nP_tilde =\n" << *P_tilde_;
	}
	M.apply(opM_notrans, *P_,       Q_.get()       );       // op(M)*P^{1}                                    -> Q^{i}
	M.apply(opM_trans,   *P_tilde_, Q_tilde_.get() );       // op(M)'*P_tilde^{1}                             -> Q_tilde^{i}
	if(get_out().get() && dump_all()) {
		*get_out() << "\nQ =\n" << *Q_;
		*get_out() << "\nQ_tilde =\n" << *Q_tilde_;
	}
	dots( *P_tilde_, *Q_, &gamma_[0] );                     // P_tilde^{i}[j]' * Q^{i}[j]                     -> gamma_{i-1}[j]
	for(j=0;j<m;++j) alpha_[j] = rho_[j]/gamma_[j];         // rho_{i-1}[j] / gamma_{i-1}[j]                  -> alpha_{i}[j]
	if(get_out().get() && dump_all()) {
		*get_out() << "\ngamma =\n"; for(j=0;j<m;++j) *get_out() << " " << gamma_[j]; *get_out() << std::endl;
		*get_out() << "\nalpha =\n"; for(j=0;j<m;++j) *get_out() << " " << alpha_[j]; *get_out() << std::endl;
	}
	for(j=0;j<m;++j) { 	// Check for failure: alpha_{i-1} = 0 or NaN or Inf
		TEST_FOR_EXCEPTION(
			alpha_[j] == 0.0 || RTOp_is_nan_inf(alpha_[j]), Exceptions::SolverBreakdown
			,"BiCGSolver<Scalar>::solve(...): Error, rho["<<j<<"] = 0.0, the method has failed!"
			);
	}
	update( &alpha_[0], +1.0, *P_, X );                    // +alpha_{i}[j] * P^{i}[j] + X^{i-1}             -> X^{i} 
	update( &alpha_[0], -1.0, *Q_, R_.get() );             // -alpha_{i}[j] * Q^{i}[j] + R^{i-1}             -> R^{i} 
	update( &alpha_[0], -1.0, *Q_tilde_, R_tilde_.get() ); // -alpha_{i}[j] * Q_tilde^{i}[j] + R_tilde^{i-1} -> R_tilde^{i} 
	if(get_out().get() && dump_all()) {
		*get_out() << "\nX =\n" << *X;
		*get_out() << "\nR =\n" << *R_;
		*get_out() << "\nR_tilde =\n" << *R_tilde_;
	}
}

template<class Scalar>
void BiCGSolver<Scalar>::compress( bool isConverged[],  MultiVector<Scalar>* X  ) const
{
	// Copy the solutions and compress out the solved systems
	Index k_curr = 0;
	for( Index k = 0; k < currNumSystems_; ++k ) {
		if( isConverged[k] ) {
			// Save the solution
			assign( X->col( activeSystems_[k] ).get(), *X_curr_->col(k+1) );
		}
		else {
			// Compress out the solved systems
			compressMember( k, k_curr, X_curr_ );
			compressMember( k, k_curr, R_ );
			compressMember( k, k_curr, R_tilde_ );
			compressMember( k, k_curr, Q_ );
			compressMember( k, k_curr, Q_tilde_ );
			compressMember( k, k_curr, Z_ );
			compressMember( k, k_curr, Z_tilde_ );
			compressMember( k, k_curr, P_ );
			compressMember( k, k_curr, P_tilde_ );
			compressMember( k, k_curr, &rho_ );
			compressMember( k, k_curr, &rho_old_ );
			compressMember( k, k_curr, &beta_ );
			compressMember( k, k_curr, &gamma_ );
			compressMember( k, k_curr, &alpha_ );
			compressMember( k, k_curr, &rel_err_denom_ );
			compressMember( k, k_curr, &norms_ );
			compressMember( k, k_curr, &activeSystems_ );
			++k_curr;
		}
	}
	currNumSystems_ = k_curr;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_BICG_SOLVER_HPP
