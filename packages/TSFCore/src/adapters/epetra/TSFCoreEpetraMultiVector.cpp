// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.cpp

#include "TSFCoreEpetraMultiVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "ThrowException.hpp"
#include "Epetra_MultiVector.h"

namespace TSFCore {

// Constructors/Initializers

EpetraMultiVector::EpetraMultiVector()
{}

EpetraMultiVector::EpetraMultiVector(
	const MemMngPack::ref_count_ptr<Epetra_MultiVector>         multi_vec
	,const MemMngPack::ref_count_ptr<const EpetraVectorSpace>   &range
	,const MemMngPack::ref_count_ptr<const EpetraVectorSpace>   &domain
	)
{
	this->initialize(multi_vec,range,domain);
}

void EpetraMultiVector::initialize(
	const MemMngPack::ref_count_ptr<Epetra_MultiVector>         multi_vec
	,const MemMngPack::ref_count_ptr<const EpetraVectorSpace>   &range
	,const MemMngPack::ref_count_ptr<const EpetraVectorSpace>   &domain
	)
{
#ifdef _DEBUG
	const char err_msg[] = "EpetraMultiVector::initialize(...): Error!";
	THROW_EXCEPTION( multi_vec.get() == NULL, std::invalid_argument, err_msg ); 
	THROW_EXCEPTION( range.get()   == NULL, std::invalid_argument, err_msg ); 
	THROW_EXCEPTION( domain.get()  == NULL, std::invalid_argument, err_msg ); 
	THROW_EXCEPTION( range->dim()  == 0,    std::invalid_argument, err_msg ); 
	THROW_EXCEPTION( domain->dim() == 0,    std::invalid_argument, err_msg );
	// ToDo: Check the compatibility of the vectors in col_vecs!
#endif
	multi_vec_ = multi_vec;
	range_     = range;
	domain_    = domain;
}

void EpetraMultiVector::setUninitialized()
{
	multi_vec_ = MemMngPack::null;
	range_     = MemMngPack::null;
	domain_    = MemMngPack::null;
}

// Overridden from LinearOp

MemMngPack::ref_count_ptr<const VectorSpace<EpetraMultiVector::Scalar> >
EpetraMultiVector::domain() const
{
	return domain_;
}

// Overridden from MultiVector

MemMngPack::ref_count_ptr<Vector<EpetraMultiVector::Scalar> >
EpetraMultiVector::col(Index j)
{
	THROW_EXCEPTION( !(  1 <= j  && j <= domain_->dim() ), std::logic_error, "EpetraMultiVector::col(j): Error!" );
	assert(0); // ToDo: Implement!
	return MemMngPack::null;
}

MemMngPack::ref_count_ptr<MultiVector<EpetraMultiVector::Scalar> >
EpetraMultiVector::subView( const Range1D& col_rng_in )
{
	return MultiVector<Scalar>::subView(col_rng_in); // ToDo: specialize!
/*
	const Index   cols    = domain_->dim();
	const Range1D col_rng = RangePack::full_range(col_rng_in,1,cols);
#ifdef _DEBUG
	THROW_EXCEPTION(
		!( col_rng.ubound() <= cols )
		,std::logic_error
		,"EpetraMultiVector::subView(col_rng): Error, the input range col_rng = ["<<col_rng.lbound()<<","<<col_rng.ubound()<<"] "
		"is not in the range [1,"<<cols<<"]!"
		);
#endif
	return MemMngPack::rcp(
		new EpetraMultiVector<Scalar>(
			range_,domain_->smallVecSpcFcty()->createVecSpc(col_rng.size()),&col_vecs_[col_rng.lbound()-1]
			) );
*/
}

// Overridden from MPIMultiVectorBase

MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<EpetraMultiVector::Scalar> >
EpetraMultiVector::mpiSpace() const
{
	return range_;
}
	
} // end namespace TSFCore
