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
// TSFCoreEpetraLinearOp.cpp

#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCore_get_Epetra_MultiVector.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

namespace TSFCore {

// Constructors / initializers / accessors

EpetraLinearOp::EpetraLinearOp()
{}

EpetraLinearOp::EpetraLinearOp(
	const Teuchos::RefCountPtr<Epetra_Operator>   &op
	,ETransp                                      opTrans
	)
{
	initialize(op,opTrans);
}

void EpetraLinearOp::initialize(
	const Teuchos::RefCountPtr<Epetra_Operator>   &op
	,ETransp                                      opTrans
	)
{
	namespace mmp = MemMngPack;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( op.get()==NULL, std::invalid_argument, "EpetraLinearOp::initialize(...): Error!" );
#endif
	op_      = op;
	opTrans_ = opTrans;
  /* KRL: call allocators to create domain and range space objects. 
   * This ensures that the correct subtype is created. */
	domain_  = allocateDomain(op, opTrans);
	range_   = allocateRange(op, opTrans);
}

void EpetraLinearOp::setUninitialized(
	Teuchos::RefCountPtr<Epetra_Operator>    *op
	,ETransp                                 *opTrans
	)
{

	if(op)      *op      = op_;
	if(opTrans) *opTrans = opTrans_;

	op_      = Teuchos::null;
	opTrans_ = NOTRANS;
	domain_  = Teuchos::null;
	range_   = Teuchos::null;

}

Teuchos::RefCountPtr<const EpetraVectorSpace>
EpetraLinearOp::epetraRange() const
{
  return range_;
}

Teuchos::RefCountPtr<const EpetraVectorSpace>
EpetraLinearOp::epetraDomain() const
{
  return domain_;
}

const Teuchos::RefCountPtr<Epetra_Operator>&
EpetraLinearOp::epetra_op() const 
{
  return op_;
}

Teuchos::RefCountPtr<Epetra_Operator>
EpetraLinearOp::epetra_op() 
{
  return op_;
}

// Overridden from OpBase

Teuchos::RefCountPtr<const VectorSpace<EpetraLinearOp::Scalar> >
EpetraLinearOp::range() const
{
	return range_;
}

Teuchos::RefCountPtr<const VectorSpace<EpetraLinearOp::Scalar> >
EpetraLinearOp::domain() const
{
	return domain_;
}


// Overridden from LinearOp

Teuchos::RefCountPtr<const LinearOp<EpetraLinearOp::Scalar> >
EpetraLinearOp::clone() const
{
	namespace mmp = MemMngPack;
	assert(0); // ToDo: Implement when needed
	return Teuchos::null;
}

void EpetraLinearOp::apply(
	const ETransp            M_trans
	,const Vector<Scalar>    &x
	,Vector<Scalar>          *y
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	MultiVectorCols<Scalar>
		X(Teuchos::rcp(const_cast<Vector<Scalar>*>(&x),false)),
		Y(Teuchos::rcp(y,false));
	apply(M_trans,X,&Y,alpha,beta);
}

void EpetraLinearOp::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X_in
	,MultiVector<Scalar>          *Y_inout
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
#ifdef _DEBUG
	// ToDo: Assert vector spaces!
#endif
	//
	// Get Epetra_MultiVector objects for the arguments
	//
	Teuchos::RefCountPtr<const Epetra_MultiVector>
		X = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *domain_ : *range_
			,Teuchos::rcp(&X_in,false)
			);
	Teuchos::RefCountPtr<Epetra_MultiVector>
		Y;
	if( beta == 0 ) {
		Y = get_Epetra_MultiVector(
			M_trans==NOTRANS ? *range_ : *domain_
			,Teuchos::rcp(Y_inout,false)
			);
	}
	//
	// Set the operator mode
	//
	op_->SetUseTranspose( trans_trans(opTrans_,M_trans) == NOTRANS ? false : true );
	//
	// Perform the operation
	//
	if( beta == 0 ) {
		// Y = M * X
		op_->Apply(
			*X
			,*Y
			);
		// Y = alpha * Y
		if( alpha != 1.0 ) Y->Scale(alpha);
	}
	else {
		// Y_inout = beta * Y_inout
		if(beta != 0.0) scale( beta, Y_inout );
		else assign( Y_inout, 0.0 );
		// T = M * X
		Epetra_MultiVector T(
			M_trans == NOTRANS ? op_->OperatorRangeMap() : op_->OperatorDomainMap()
			,X_in.domain()->dim()
			,false
			);
		op_->Apply(
			*X
			,T
			);
		// Y_inout += alpha * T
		update(
      alpha
      ,EpetraMultiVector(
        Teuchos::rcp( &T, false)
        ,Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(Y_inout->range())
        ,Teuchos::rcp_dynamic_cast<const EpetraVectorSpace>(Y_inout->domain())
        )
      ,Y_inout
      );
	}
}


// --------- KRL: allocator methods for domain and range spaces ---- 

Teuchos::RefCountPtr<const EpetraVectorSpace> 
EpetraLinearOp::allocateDomain(const Teuchos::RefCountPtr<Epetra_Operator>  &op 
                               ,ETransp  op_trans 
                               )  const
{
  return Teuchos::rcp( new EpetraVectorSpace( Teuchos::rcp(&op->OperatorDomainMap(),false) ) );
}

Teuchos::RefCountPtr<const EpetraVectorSpace> 
EpetraLinearOp::allocateRange(const Teuchos::RefCountPtr<Epetra_Operator>  &op 
                              ,ETransp  op_trans 
                               )  const
{
  return Teuchos::rcp( new EpetraVectorSpace( Teuchos::rcp(&op->OperatorRangeMap(),false) ) );
}



}	// end namespace TSFCore
