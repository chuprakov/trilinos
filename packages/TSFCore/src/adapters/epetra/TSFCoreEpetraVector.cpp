// /////////////////////////////////////////////////////////////////
// TSFCoreEpetraVector.cpp

#include "TSFCoreEpetraVector.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Teuchos_TestForException.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace TSFCore {

// Constructors/initializers

EpetraVector::EpetraVector()
{}

EpetraVector::EpetraVector(
	const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
	)
{
	initialize(epetra_vec,epetra_vec_spc);
}

void EpetraVector::initialize(
	const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
	,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
	)
{
  using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_vec.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
	TEST_FOR_EXCEPTION( !epetra_vec_spc.get(), std::invalid_argument, "EpetraVector::initialize(...): Error!");
#endif
	epetra_vec_ = epetra_vec;
  epetra_vec_spc_ = epetra_vec_spc;
	updateMpiSpace();
}

void EpetraVector::setUninitialized(
	Teuchos::RefCountPtr<Epetra_Vector>              *epetra_vec
	,Teuchos::RefCountPtr<const EpetraVectorSpace>   *epetra_vec_spc
	)
{
	if(epetra_vec) *epetra_vec = epetra_vec_;
	if(epetra_vec_spc) *epetra_vec_spc = epetra_vec_spc_;

	epetra_vec_ = Teuchos::null;
	epetra_vec_spc_ = Teuchos::null;

	updateMpiSpace();
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
