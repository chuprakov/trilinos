// /////////////////////////////////////////////////////////////////
// TSFCoreEpetraVector.cpp

#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

#include "Epetra_Vector.h"

namespace TSFCore {

// Constructors/initializers

EpetraVector::EpetraVector()
{}

EpetraVector::EpetraVector(
	const Teuchos::RefCountPtr<Epetra_Vector>             &epetra_vec
	)
{
	initialize(epetra_vec);
}

void EpetraVector::initialize(
	const Teuchos::RefCountPtr<Epetra_Vector>             &epetra_vec
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_vec.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
#endif
	//
	epetra_vec_     = epetra_vec;
	epetra_vec_spc_ = Teuchos::rcp( new EpetraVectorSpace( Teuchos::rcp( &epetra_vec->Map(), false ) ) );
}

Teuchos::RefCountPtr<Epetra_Vector>
EpetraVector::setUninitialized()
{

	Teuchos::RefCountPtr<Epetra_Vector> tmp = epetra_vec_;
	
	epetra_vec_ = Teuchos::null;
	epetra_vec_spc_ = Teuchos::null;

	return tmp;
}

// Overridden from MPIVectorBase

Teuchos::RefCountPtr<const MPIVectorSpaceBase<EpetraVector::Scalar> >
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
