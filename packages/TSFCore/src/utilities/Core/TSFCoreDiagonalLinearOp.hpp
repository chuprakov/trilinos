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
// TSFCoreDiagonalLinearOp.hpp

#ifndef TSFCORE_DIAGONAL_LINEAR_OP_HPP
#define TSFCORE_DIAGONAL_LINEAR_OP_HPP

#include "TSFCoreDiagonalLinearOpDecl.hpp"
#include "TSFCoreVector.hpp"

namespace TSFCore {

// Constructors/initializers/accessors

template<class Scalar>
DiagonalLinearOp<Scalar>::DiagonalLinearOp()
{}

template<class Scalar>
DiagonalLinearOp<Scalar>::DiagonalLinearOp(
	const Teuchos::RefCountPtr<const Vector<Scalar> >   &diag
	)
{
	initialize(diag);
}

template<class Scalar>
void DiagonalLinearOp<Scalar>::initialize(
	const Teuchos::RefCountPtr<const Vector<Scalar> >   &diag
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPT(diag.get()==NULL);
#endif
	diag_ = diag;
}

template<class Scalar>
void DiagonalLinearOp<Scalar>::uninitialize(
	Teuchos::RefCountPtr<const Vector<Scalar> >  *diag
	)
{
	if(diag) *diag = diag_;
	diag_ = Teuchos::null;
}

// Overridden from OpBase

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
DiagonalLinearOp<Scalar>::domain() const
{
	return diag_->space();
}

template<class Scalar>
Teuchos::RefCountPtr< const VectorSpace<Scalar> >
DiagonalLinearOp<Scalar>::range() const
{
	return diag_->space();
}

template<class Scalar>
bool DiagonalLinearOp<Scalar>::opSupported(ETransp M_trans) const
{
	return true;
}

// Overridden from LinearOp

template<class Scalar>
void DiagonalLinearOp<Scalar>::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	typedef Teuchos::ScalarTraits<Scalar> ST;
	if( beta != ST::one() ) Vt_S( y, beta );
	ele_wise_prod( alpha, x, *diag_, y );
}

template<class Scalar>
void DiagonalLinearOp<Scalar>::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	LinearOp<Scalar>::apply(M_trans,X,Y,alpha,beta); // ToDo: specialize if we can?
}

template<class Scalar>
Teuchos::RefCountPtr<const LinearOp<Scalar> >
DiagonalLinearOp<Scalar>::clone() const
{
	return Teuchos::null; // Not supported yet but could be
}

}	// end namespace TSFCore

#endif	// TSFCORE_DIAGONAL_LINEAR_OP_HPP
