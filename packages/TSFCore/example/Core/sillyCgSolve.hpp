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

#ifndef TSFCORE_SILLY_CG_SOLVE_HPP
#define TSFCORE_SILLY_CG_SOLVE_HPP

#include "TSFCoreLinOp.hpp"
#include "TSFCoreVectorStdOps.hpp"

///
/** \brief Silly little example unpreconditioned CG solver
 *
 * This little function is just a silly little ANA that implements the
 * CG (conjugate gradient) method for solving symmetric positive definite
 * systems using the \ref fundamental_TSFCore_interfaces_sec "fundamental TSFCore interfaces".
 *
 * This function is small and is ment to be looked at so study its
 * implementation by clicking on the below link to its definition.
 *
 * \ingroup TSFCore_examples_cg_grp
 */
template<class Scalar>
bool sillyCgSolve(
	const TSFCore::LinOpNonPersisting<Scalar>                      &A
	,const TSFCore::Vector<Scalar>                                 &b
	,const int                                                     maxNumIters
	,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType   tolerance
	,TSFCore::Vector<Scalar>                                       *x
	,std::ostream                                                  *out          = NULL
	)
{
	using Teuchos::RefCountPtr;
	typedef Teuchos::ScalarTraits<Scalar>   ST;           // We need to use ScalarTraits to support arbitrary types!         
	typedef typename ST::magnitudeType      ScalarMag;    // This is the type for a norm
	// Validate input
	TEST_FOR_EXCEPT(x==NULL);
	TSFCORE_ASSERT_LINEAR_OP_VEC_APPLY_SPACES("sillyCgSolve()",*A.op(),A.defaultTrans(),*x,&b); // A*x - b agree?
	if(out) *out << "\nStarting CG solver ...\n" << std::scientific << "\n  Type of A = \'"<<typeid(*A.op()).name()<<"\'"
							 << "\n  Type of b = \'"<<typeid(b).name()<<"\'" << "\n  Type of x = \'"<<typeid(*x).name()<<"\'\n\n";
	// Get the vector space (domain and range spaces should be the same)
	RefCountPtr<const TSFCore::VectorSpace<Scalar> > space = A.domain();
	// Compute initial residual : r = b - A*x
	RefCountPtr<TSFCore::Vector<Scalar> > r = space->createMember();                     
	TSFCore::assign(&*r,b);                               // r = b
	A.apply(TSFCore::NOTRANS,*x,&*r,-ST::one(),ST::one());// r = -A*x + r
	const ScalarMag r0_nrm = TSFCore::norm(*r);           // Compute ||r0|| = sqrt(<r0,r0>) for convergence test
	if(r0_nrm == ST::zero()) return true;                 // Trivial RHS and initial LHS guess?
	// Create workspace vectors and scalars
	RefCountPtr<TSFCore::Vector<Scalar> > p = space->createMember(), q = space->createMember();
	Scalar rho_old;
	// Perform the iterations
	for( int iter = 0; iter <= maxNumIters; ++iter ) {
		// Check convergence and output iteration
		const ScalarMag r_nrm = TSFCore::norm(*r);          // Compute ||r|| = sqrt(<r,r>)
		const bool isConverged = r_nrm/r0_nrm <= tolerance;
		if( iter%(maxNumIters/10+1) == 0 || iter == maxNumIters || isConverged ) {
			if(out) *out << "Iter = " << iter << ", ||b-A*x||/||b-A*x0|| = " << (r_nrm/r0_nrm) << std::endl;
			if( r_nrm/r0_nrm < tolerance ) return true;       // Converged to tolerance, Success!
		}
		// Compute iteration
		const Scalar rho = space->scalarProd(*r,*r);        // <r,r>              -> rho
		if(iter==0) TSFCore::assign(&*p,*r);                // r                  -> p   (iter == 0)
		else TSFCore::Vp_V( &*p, *r, Scalar(rho/rho_old) ); // r+(rho/rho_old)*p  -> p   (iter  > 0)
		A.apply(TSFCore::NOTRANS,*p,&*q);                   // A*p                -> q
		const Scalar alpha = rho/space->scalarProd(*p,*q);  // rho/<p,q>          -> alpha
		TSFCore::Vp_StV( x,   Scalar(+alpha), *p );         // +alpha*p + x       -> x
		TSFCore::Vp_StV( &*r, Scalar(-alpha), *q );         // -alpha*q + r       -> r
		rho_old = rho;                                      // rho                -> rho_old (remember rho for next iter)
	}
	return false; // Failure
} // end sillyCgSolve

#endif // TSFCORE_SILLY_CG_SOLVE_HPP
