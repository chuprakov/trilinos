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

// //////////////////////////////////////////////////////////////////////////
// TSFCoreSolversGmresSolver.hpp
//

#ifndef TSFCORE_SOLVERS_GMRES_SOLVER_HPP
#define TSFCORE_SOLVERS_GMRES_SOLVER_HPP

#include "TSFCoreSolversTypes.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "TSFCoreMultiVector.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_BLAS.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
class GMRESSolver {

public:

	GMRESSolver(
		int		    default_max_iter = 1000,
		Scalar		default_tol      = 1e-10			
		);

	int currIteration() const { return curr_iter; };
  
	Scalar currEstRelResidualNorm() const { return curr_res; };
  
	SolveReturn solve(
		const LinearOp<Scalar> &Op,
		const Vector<Scalar>   &b,
		Vector<Scalar>         *curr_soln,
		const ETransp          trans_in,
		int                    max_iter_in,
		Scalar                 tol_in 
		);	

	Teuchos::RefCountPtr<const GMRESSolver<Scalar> > clone () const;
	
private:

	void doIteration( const LinearOp<Scalar> &Op, const ETransp Op_trans );
  
	int				                                              max_iter, curr_iter;
	bool				                                          isConverged;
	Scalar			                                              tol, curr_res, r0; 
	std::vector<Scalar>                                           z;
	Teuchos::SerialDenseMatrix<int,Scalar>                        H_;
	Teuchos::RefCountPtr< MultiVector<Scalar> >                   V_;
	Teuchos::RefCountPtr< Vector<Scalar> >                        r;
	std::vector<Scalar>		                                      cs, sn;	

}; // end class GMRESSolver

template<class Scalar>
GMRESSolver<Scalar>::GMRESSolver(
	int			default_max_iter,
	Scalar		default_tol
	)
    : max_iter(default_max_iter), tol(default_tol),
      isConverged(false), curr_iter(0), r0( 1.0 ), curr_res( 1.0 )
{}
	
template<class Scalar>
SolveReturn GMRESSolver<Scalar>::solve(
	const LinearOp<Scalar> &Op,
	const Vector<Scalar>   &b,
	Vector<Scalar>         *curr_soln,
	const ETransp          Op_trans,
	int                    max_iter_in,
	Scalar                 tol_in
	)
{
    // 
    //  Check compatability of linear operator with rhs and solution vector.
    //
    const VectorSpace<Scalar> &Op_domain = ( ( Op_trans == NOTRANS ) ? *Op.domain() : *Op.range() );
    const VectorSpace<Scalar> &Op_range = ( ( Op_trans == NOTRANS ) ? *Op.range() : *Op.domain() );
    const bool 	domain_compatable = Op_domain.isCompatible( *curr_soln->space() ),
		range_compatable = Op_range.isCompatible( *b.space() );
    if (!(domain_compatable && range_compatable) ) { std::cout<<"Op is not compatable with x or b"<<std::endl; }
    //
    //  Initialize internal data structures
    //		
	curr_iter = 0;
    V_ = Op_domain.createMembers(max_iter+1);
    r = Op_domain.createMember();
    //
    if ( tol_in != tol ) { tol = tol_in; }
	if( max_iter_in+1 != H_.numRows() ) {
		max_iter = max_iter_in;
		H_.shape( max_iter+1, max_iter );
		z.resize(max_iter+1);
		cs.resize(max_iter); sn.resize(max_iter);
    }
	//
	// Check if the RHS is zero
	//
	isConverged = false;
	const Scalar norm_2_b = norm_2( b );
    if(norm_2_b == 0.0) {
		isConverged = true;
		assign( curr_soln, 0.0 );
	}
	if(!isConverged) {
		//
		// Determine the residual from the current solution: r = b - Op*curr_soln 
		//
		assign( r.get(), b );
		Op.apply( Op_trans, *curr_soln, r.get(), 1.0, -1.0 );
		curr_res = norm_2( *r ) / ( 1.0 + norm_2( b ) );
		if (curr_res < tol) isConverged = true;
		//
		// Set up initial vector.
		//
		r0 = norm_2( *r );
		z[0] = r0;
		Teuchos::RefCountPtr<Vector<Scalar> >
			w = V_->col(1);		    // get a mutable view of the first column of V_.
		assign( w.get(), *r );      // copy r to the first column of V_.
		Vt_S( w.get(), 1.0/r0 );    // v_1 = r_0 / ||r_0||
		w = Teuchos::null;       // Forces V_ to be modified
		//
		// Calls doIteration() using the current linear system parameters.
		//
		while( !isConverged && (curr_iter < max_iter) ) { doIteration( Op, Op_trans ); }
		//
		// Solve least squares problem.
		//
		Teuchos::BLAS<int, Scalar> blas;
		blas.TRSM(Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS, Teuchos::NON_UNIT_DIAG, 
				  curr_iter, 1, 1.0, H_.values(), max_iter+1, &z[0], max_iter+1 );
		//
		// Compute the new solution.
		//
		if( curr_iter > 0 ) {
			Teuchos::RefCountPtr<MultiVector<Scalar> > V = V_->subView(Range1D(1,curr_iter));
			SerialVector<Scalar> z_vec( &z[0], 1, curr_iter, false );
			V->apply( NOTRANS, z_vec, curr_soln, -1.0, 1.0 );
		}
		//for( int i = 0; i < curr_iter; i++ )
		//	Vp_StV( curr_soln, -z[i], *V_->col(i+1) );
	}
	return SolveReturn( curr_iter >= max_iter ? MAX_ITER_EXCEEDED : SOLVED_TO_TOL, curr_iter );
}

template<class Scalar>
Teuchos::RefCountPtr<const GMRESSolver<Scalar> >
GMRESSolver<Scalar>::clone() const
{
	return Teuchos::rcp( new GMRESSolver<Scalar>(*this) );
}

template<class Scalar>
void GMRESSolver<Scalar>::doIteration( const LinearOp<Scalar> &Op, const ETransp Op_trans )
{
    int i;
    Scalar temp;
    Teuchos::BLAS<int, Scalar> blas;
    Teuchos::SerialDenseMatrix<int, Scalar> &H = H_;
    // 
	Teuchos::RefCountPtr<Vector<Scalar> >
		w = V_->col(curr_iter+2);                              // w = v_{j+1}
    Op.apply( Op_trans, *V_->col(curr_iter+1), w.get() );      // w = Op * v_{j}	
    //
    // Perform MGS to orthogonalize new Krylov vector.
    //
   for( i=0; i<curr_iter+1; i++ ) {	
		H( i, curr_iter ) = dot( *w, *V_->col(i+1) );          // h_{i,j} = ( w, v_{i} )
		Vp_StV( w.get(), -H( i, curr_iter ), *V_->col(i+1) );  // w = w - h_{i,j} * v_{i}
    }
/* RAB: Why not this or a MultiVector version???
    for( i=0; i<curr_iter+1; i++ ) {
		H( i, curr_iter ) = dot( *w, *V_->col(i+1) );          // h_{i,j} = ( w, v_{i} )
    }
    for( i=0; i<curr_iter+1; i++ ) {
		Vp_StV( w.get(), -H( i, curr_iter ), *V_->col(i+1) );  // w = w - h_{i,j} * v_{i}
    }
*/
    H( curr_iter+1, curr_iter ) = norm_2( *w );                // h_{j+1,j} = || w ||
    Vt_S( w.get(), 1.0 / H( curr_iter+1, curr_iter ) ); 	   // v_{j+1} = w / h_{j+1,j}			
    //
    // Apply previous Givens rotations
    //
    for( i=0; i<curr_iter; i++ ) {
		temp = cs[i]*H( i, curr_iter ) + sn[i]*H( i+1, curr_iter );
		H( i+1, curr_iter ) = -sn[i]*H( i, curr_iter ) + cs[i]*H( i+1, curr_iter );
		H( i, curr_iter ) = temp;
    }
    //
    // Calculate new Givens rotation
    //
    blas.ROTG( &H( curr_iter, curr_iter ), &H( curr_iter+1, curr_iter ), 
			   &cs[curr_iter], &sn[curr_iter] );
    //
    // Update RHS and residual w/ new transform and compute residual.
    //
    z[curr_iter+1] = -sn[curr_iter]*z[curr_iter];
    z[curr_iter] *= cs[curr_iter];
    curr_res = Teuchos::ScalarTraits<Scalar>::magnitude( z[curr_iter+1] ) / r0; 
    if (curr_res < tol) { isConverged = true; }
    //    
    // Increment the iteration counter.
    //
    curr_iter++;
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_GMRES_SOLVER_HPP
