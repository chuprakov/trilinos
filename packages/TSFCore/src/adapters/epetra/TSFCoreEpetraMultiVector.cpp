// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.cpp

#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_MultiVector.h"

namespace TSFCore {

// Constructors/Initializers

EpetraMultiVector::EpetraMultiVector()
{}

EpetraMultiVector::EpetraMultiVector(
	const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain
	)
{
	this->initialize(epetra_multi_vec,epetra_range,epetra_domain);
}

void EpetraMultiVector::initialize(
	const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain
	)
{
#ifdef _DEBUG
	const char err_msg[] = "EpetraMultiVector::initialize(...): Error!";
	TEST_FOR_EXCEPTION( epetra_multi_vec.get() == NULL, std::invalid_argument, err_msg ); 
	TEST_FOR_EXCEPTION( epetra_range.get() && epetra_range->dim() == 0, std::invalid_argument, err_msg ); 
	TEST_FOR_EXCEPTION( epetra_domain.get() && epetra_domain->dim() == 0,std::invalid_argument, err_msg );
	// ToDo: Check the compatibility of the vectors in col_vecs!
#endif
	epetra_multi_vec_ = epetra_multi_vec;
	if(epetra_range.get()) {
		epetra_range_  = epetra_range;
	}
	else {
		assert(0); // ToDo: Implement this case!
	}
	if(epetra_domain.get()) {
		epetra_domain_  = epetra_domain;
	}
	else {
		assert(0); // ToDo: Implement this case!
	}
}

void EpetraMultiVector::setUninitialized(
	Teuchos::RefCountPtr<Epetra_MultiVector>        *epetra_multi_vec
	,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_range
	,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_domain
	)
{
	if(epetra_multi_vec) *epetra_multi_vec = epetra_multi_vec_;;
	if(epetra_range) *epetra_range = epetra_range_;
	if(epetra_domain) * epetra_domain = epetra_domain_;
	epetra_multi_vec_ = Teuchos::null;
	epetra_range_ = Teuchos::null;
	epetra_domain_ = Teuchos::null;
}

// Overridden from OpBase

Teuchos::RefCountPtr<const VectorSpace<EpetraMultiVector::Scalar> >
EpetraMultiVector::domain() const
{
	return epetra_domain_;
}

// Overridden from LinearOp

void EpetraMultiVector::apply(
	const ETransp                 M_trans
	,const MultiVector<Scalar>    &X
	,MultiVector<Scalar>          *Y
	,const Scalar                 alpha
	,const Scalar                 beta
	) const
{
	assert(0); // ToDo: Implement this!
}

// Overridden from MultiVector

Teuchos::RefCountPtr<Vector<EpetraMultiVector::Scalar> >
EpetraMultiVector::col(Index j)
{
	TEST_FOR_EXCEPTION( !(  1 <= j  && j <= epetra_domain_->dim() ), std::logic_error, "EpetraMultiVector::col(j): Error!" );
	assert(0); // ToDo: Implement!
	return Teuchos::null;
}

Teuchos::RefCountPtr<MultiVector<EpetraMultiVector::Scalar> >
EpetraMultiVector::subView( const Range1D& col_rng_in )
{
	return MultiVector<Scalar>::subView(col_rng_in); // ToDo: specialize!
/*
	const Index   cols    = domain_->dim();
	const Range1D col_rng = RangePack::full_range(col_rng_in,1,cols);
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		!( col_rng.ubound() <= cols )
		,std::logic_error
		,"EpetraMultiVector::subView(col_rng): Error, the input range col_rng = ["<<col_rng.lbound()<<","<<col_rng.ubound()<<"] "
		"is not in the range [1,"<<cols<<"]!"
		);
#endif
	return Teuchos::rcp(
		new EpetraMultiVector<Scalar>(
			range_,domain_->smallVecSpcFcty()->createVecSpc(col_rng.size()),&col_vecs_[col_rng.lbound()-1]
			) );
*/
}

Teuchos::RefCountPtr<MultiVector<EpetraMultiVector::Scalar> >
EpetraMultiVector::subView( const int numCols, const int cols[] )
{
	return MultiVector<Scalar>::subView(numCols,cols); // ToDo: specialize!
}

// Overridden from MPIMultiVectorBase

Teuchos::RefCountPtr<const MPIVectorSpaceBase<EpetraMultiVector::Scalar> >
EpetraMultiVector::mpiSpace() const
{
	return epetra_range_;
}

void EpetraMultiVector::getLocalData( Scalar **values, Index *leadingDim )
{
	assert(0); // ToDo: Implement!
}
	
} // end namespace TSFCore
