// /////////////////////////////////////////////////////////////////
// TSFCoreEpetraVector.cpp

#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "ThrowException.hpp"

#include "Epetra_Vector.h"

namespace TSFCore {

// Constructors/initializers

EpetraVector::EpetraVector()
	:globalDim_(-1),localOffset_(-1),localSubDim_(-1)
{}

EpetraVector::EpetraVector(
	const MemMngPack::ref_count_ptr<Epetra_Vector>             &epetra_vec
	)
{
	initialize(epetra_vec);
}

void EpetraVector::initialize(
	const MemMngPack::ref_count_ptr<Epetra_Vector>             &epetra_vec
	)
{
	namespace mmp = MemMngPack;
#ifdef _DEBUG
	THROW_EXCEPTION( !epetra_vec.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
#endif
	//
	epetra_vec_     = epetra_vec;
	epetra_vec_spc_ = mmp::rcp( new EpetraVectorSpace( mmp::rcp( &epetra_vec->Map(), false ) ) );
	// Cache some stuff
	globalDim_    = epetra_vec_spc_->dim();
	localOffset_  = epetra_vec_spc_->localOffset();
	localSubDim_  = epetra_vec_spc_->localSubDim();
}

MemMngPack::ref_count_ptr<Epetra_Vector>
EpetraVector::setUninitialized()
{
	namespace mmp = MemMngPack;

	mmp::ref_count_ptr<Epetra_Vector> tmp = epetra_vec_;
	
	epetra_vec_ = mmp::null;
	epetra_vec_spc_ = mmp::null;
	globalDim_      = -1;
	localOffset_    = -1;
	localSubDim_    = -1;

	return tmp;
}
	
// Overridden from Vector

void EpetraVector::getSubVector( const Range1D& rng_in, RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	const Range1D rng = validateRange(rng_in);
	if( rng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
		return;
	}
	// rng consists of all local data so get it!
	const Scalar *localValues = &(*epetra_vec_)[0];
	sub_vec->initialize(
		rng.lbound()-1                             // globalOffset
		,rng.size()                                // subDim
		,localValues+(rng.lbound()-localOffset_-1) // values
		,1                                         // stride
		);
}

void EpetraVector::freeSubVector( RTOpPack::SubVectorT<Scalar>* sub_vec ) const
{
	if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
		// Let the default implementation handle it!
		Vector<Scalar>::freeSubVector(sub_vec);
		return;
	}
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

void EpetraVector::getSubVector( const Range1D& rng_in, RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	const Range1D rng = validateRange(rng_in);
	if( rng.lbound() < localOffset_+1 || localOffset_+localSubDim_ < rng.ubound() ) {
		// rng consists of off-processor elements so use the default implementation!
		Vector<Scalar>::getSubVector(rng_in,sub_vec);
		return;
	}
	// rng consists of all local data so get it!
	Scalar *localValues = &(*epetra_vec_)[0];
	sub_vec->initialize(
		rng.lbound()-1                             // globalOffset
		,rng.size()                                // subDim
		,localValues+(rng.lbound()-localOffset_-1) // values
		,1                                         // stride
		);
}

void EpetraVector::commitSubVector( RTOpPack::MutableSubVectorT<Scalar>* sub_vec )
{
	if( sub_vec->globalOffset() < localOffset_ || localOffset_+localSubDim_ < sub_vec->globalOffset()+sub_vec->subDim() ) {
		// Let the default implementation handle it!
		Vector<Scalar>::commitSubVector(sub_vec);
		return;
	}
	sub_vec->set_uninitialized();  // Nothing to deallocate!
}

// Overridden from MPIVectorBase

MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<EpetraVector::Scalar> >
EpetraVector::mpiSpace() const
{
	return epetra_vec_spc_; // will be null if this is uninitialized!
}

// private

Range1D EpetraVector::validateRange( const Range1D &rng_in ) const
{
	const Range1D rng = RangePack::full_range(rng_in,1,globalDim_);
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.lbound() < 1 || globalDim_ < rng.ubound(), std::invalid_argument
		,"EpetraVector::getLocalData(...): Error, the range ["<<rng.lbound()<<","<<rng.ubound()<<"] is not "
		"in the range [1,"<<globalDim_<<"]!"
		);
#endif
	return rng;
}

} // end namespace TSFCore
