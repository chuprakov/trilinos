// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraLinearOp.cpp

#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "dynamic_cast_verbose.hpp"
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
	domain_  = Teuchos::rcp( new EpetraVectorSpace( Teuchos::rcp(&op->OperatorDomainMap(),false) ) );
	range_   = Teuchos::rcp( new EpetraVectorSpace( Teuchos::rcp(&op->OperatorRangeMap(),false) ) );
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
	,const Vector<Scalar>    &x_in
	,Vector<Scalar>          *y_inout
	,const Scalar            alpha
	,const Scalar            beta
	) const
{
	namespace mmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;

	Vp_MtV_assert_compatibility(y_inout,*this,M_trans,x_in);

	const EpetraVector   &x = dyn_cast<const EpetraVector>(x_in);
	EpetraVector         &y = dyn_cast<EpetraVector>(*y_inout);
	
	op_->SetUseTranspose( trans_trans(opTrans_,M_trans) == NOTRANS ? false : true );

	if( alpha == 1.0 && beta == 0 ) {
		Teuchos::RefCountPtr<Epetra_Vector> epetra_y = y.setUninitialized();
		op_->Apply(
			*x.epetra_vec()
			,*epetra_y
			);
		y.initialize(epetra_y);
	}
	else {
		Vt_S( y_inout, beta );
		Epetra_Vector t( M_trans == NOTRANS ? op_->OperatorRangeMap() : op_->OperatorDomainMap() );
		op_->Apply(
			*x.epetra_vec()
			,t
			);
		Vp_StV( y_inout, alpha, EpetraVector( Teuchos::rcp( &t, false) ) );
	}

}

}	// end namespace TSFCore
