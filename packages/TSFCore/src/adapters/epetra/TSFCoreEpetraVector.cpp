// /////////////////////////////////////////////////////////////////
// TSFCoreEpetraVector.cpp

#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "ThrowException.hpp"

#include "Epetra_Vector.h"

namespace TSFCore {

// Constructors/initializers

EpetraVector::EpetraVector()
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
}

MemMngPack::ref_count_ptr<Epetra_Vector>
EpetraVector::setUninitialized()
{
	namespace mmp = MemMngPack;

	mmp::ref_count_ptr<Epetra_Vector> tmp = epetra_vec_;
	
	epetra_vec_ = mmp::null;
	epetra_vec_spc_ = mmp::null;

	return tmp;
}

// Overridden from MPIVectorBase

MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<EpetraVector::Scalar> >
EpetraVector::mpiSpace() const
{
	return epetra_vec_spc_; // will be null if this is uninitialized!
}

void EpetraVector::getLocalData( Scalar** values, ptrdiff_t* stride )
{
	*values = &(*epetra_vec_)[0];
	*stride = 1;
}

} // end namespace TSFCore
