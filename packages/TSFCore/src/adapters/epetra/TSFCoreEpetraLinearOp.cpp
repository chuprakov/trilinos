// //////////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraLinearOp.cpp

#include "TSFCoreEpetraLinearOp.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"

namespace TSFCore {

// Constructors / initializers / accessors

EpetraLinearOp::EpetraLinearOp()
{}

EpetraLinearOp::EpetraLinearOp(
	const MemMngPack::ref_count_ptr<Epetra_Operator>                &op
	,ETransp                                                        opTrans
	,const MemMngPack::ref_count_ptr<MemMngPack::ReleaseResource>   &extra
	)
{
	initialize(op,opTrans,extra);
}

void EpetraLinearOp::initialize(
	const MemMngPack::ref_count_ptr<Epetra_Operator>                &op
	,ETransp                                                        opTrans
	,const MemMngPack::ref_count_ptr<MemMngPack::ReleaseResource>   &extra
	)
{
	namespace mmp = MemMngPack;
#ifdef _DEBUG
	THROW_EXCEPTION( op.get()==NULL, std::invalid_argument, "EpetraLinearOp::initialize(...): Error!" );
#endif
	state_.op      = op;
	state_.opTrans = opTrans;
	state_.extra   = extra;
	domain_        = mmp::rcp( new EpetraVectorSpace( mmp::rcp(&op->OperatorDomainMap(),false) ) );
	range_         = mmp::rcp( new EpetraVectorSpace( mmp::rcp(&op->OperatorRangeMap(),false) ) );
}

EpetraLinearOp::State
EpetraLinearOp::setUninitialized()
{
	namespace mmp = MemMngPack;

	State tmp = state_;
	
	state_   = State();
	domain_  = mmp::null;
	range_   = mmp::null;

	return tmp;
}

// Overridden from OpBase

MemMngPack::ref_count_ptr<const VectorSpace<EpetraLinearOp::Scalar> >
EpetraLinearOp::range() const
{
	return range_;
}

MemMngPack::ref_count_ptr<const VectorSpace<EpetraLinearOp::Scalar> >
EpetraLinearOp::domain() const
{
	return domain_;
}

// Overridden from LinearOp

MemMngPack::ref_count_ptr<const LinearOp<EpetraLinearOp::Scalar> >
EpetraLinearOp::clone() const
{
	namespace mmp = MemMngPack;
	assert(0);
	return mmp::null;
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
	
	state_.op->SetUseTranspose( trans_trans(state_.opTrans,M_trans) == NOTRANS ? false : true );

	if( alpha == 1.0 && beta == 0 ) {
		mmp::ref_count_ptr<Epetra_Vector> epetra_y = y.setUninitialized();
		state_.op->Apply(
			*x.epetra_vec()
			,*epetra_y
			);
		y.initialize(epetra_y);
	}
	else {
		Vt_S( y_inout, beta );
		Epetra_Vector t( M_trans == NOTRANS ? state_.op->OperatorRangeMap() : state_.op->OperatorDomainMap() );
		state_.op->Apply(
			*x.epetra_vec()
			,t
			);
		Vp_StV( y_inout, alpha, EpetraVector( mmp::rcp( &t, false) ) );
	}

}

}	// end namespace TSFCore
