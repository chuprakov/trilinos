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
// TSFCoreMultiplicativeLinearOp.hpp

#ifndef TSFCORE_MULTIPLICATIVE_LINEAR_OP_HPP
#define TSFCORE_MULTIPLICATIVE_LINEAR_OP_HPP

#include "TSFCoreMultiplicativeLinearOpDecl.hpp"
#include "TSFCoreAssertOp.hpp"

namespace TSFCore {

// Constructors/initializers/accessors

template<class Scalar>
MultiplicativeLinearOp<Scalar>::MultiplicativeLinearOp()
{}

template<class Scalar>
MultiplicativeLinearOp<Scalar>::MultiplicativeLinearOp(
	const int                        numOps
	,const LinOpPersisting<Scalar>   Ops[]
	,const Scalar                    &gamma
	)
{
	initialize(numOps,Ops,gamma);
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::initialize(
	const int                        numOps
	,const LinOpPersisting<Scalar>   Ops[]
	,const Scalar                    &gamma
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( numOps <= 0 || Ops == NULL );
	for( int k = 0; k < numOps; ++k ) {
		TEST_FOR_EXCEPT( Ops[k].op().get() == NULL );
		if( k < numOps-1 ) {
			TSFCORE_ASSERT_VEC_SPACES(
				"MultiplicativeLinearOp<Scalar>::initialize(...)"
				,*Ops[k].domain(), *Ops[k+1].range()
				);
		}
	}
#endif
	Ops_.resize(numOps);
	std::copy( Ops, Ops+numOps, Ops_.begin() );
	gamma_ = gamma;
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::uninitialize(
	const int                  numOps
	,LinOpPersisting<Scalar>   Ops[]
	,Scalar                    *gamma
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT( Ops!=NULL && numOps!=this->numOps() );
#endif
	if(Ops) std::copy( Ops_.begin(), Ops_.end(), Ops );
	if(gamma) *gamma = gamma_;
	Ops_.resize(0);
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MultiplicativeLinearOp<Scalar>::domain() const
{
	assertInitialized();
	return Ops_[numOps()-1].domain();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
MultiplicativeLinearOp<Scalar>::range() const
{
	assertInitialized();
	return Ops_[0].range();
}

template<class Scalar>
bool MultiplicativeLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
	bool opSupported = true;
	for( int k = 0; k < static_cast<int>(Ops_.size()); ++k )
		if(!Ops_[k].opSupported(M_trans)) opSupported = false;
	return opSupported;
	// ToDo: Cache these?
}

// Overridden from LinearOp

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	this->apply(
		M_trans
		,static_cast<const MultiVector<Scalar>&>(x)
		,static_cast<MultiVector<Scalar>*>(y)
		,alpha
		,beta
		);
}

template<class Scalar>
void MultiplicativeLinearOp<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	using Teuchos::RefCountPtr;
	using Teuchos::rcp;
	const int numOps = Ops_.size();
	const Index m = X.domain()->dim();
	if( M_trans == NOTRANS ) {
		//
		// Y = alpha * M * X + beta*Y
		// =>
		// Y = (alpha*gamma) * op(Op[0]) * op(Op[1]) * ... * op(Op[numOps-1]) * X + beta*Y
		//
		RefCountPtr<MultiVector<Scalar> > T_kp1, T_k; // Temporary propogated between loops 
		for( int k = numOps-1; k >= 0; --k ) {
			RefCountPtr<MultiVector<Scalar> >         Y_k;
			RefCountPtr<const MultiVector<Scalar> >   X_k;
			if(k==0)        Y_k = rcp(Y,false);  else Y_k = T_k = Ops_[k].range()->createMembers(m);
			if(k==numOps-1) X_k = rcp(&X,false); else X_k = T_kp1;
			if( k > 0 )
				Ops_[k].apply(NOTRANS,*X_k,&*Y_k);
			else
				Ops_[0].apply(NOTRANS,*X_k,&*Y_k,(alpha*gamma_),beta);
			T_kp1 = T_k;
		}
	}
	else {
		TEST_FOR_EXCEPT(true); // ToDo: Implement!
	}
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
MultiplicativeLinearOp<Scalar>::clone() const
{
	return Teuchos::null; // Not supported yet but could be
}

}	// end namespace TSFCore

#endif	// TSFCORE_MULTIPLICATIVE_LINEAR_OP_HPP
