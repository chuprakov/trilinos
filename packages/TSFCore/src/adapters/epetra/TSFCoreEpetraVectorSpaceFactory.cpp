// ///////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpaceFactory.cpp

#include "TSFCoreEpetraVectorSpaceFactory.hpp"
#include "TSFCoreEpetraVectorSpace.hpp"
#include "Epetra_LocalMap.h"
#include "Epetra_Comm.h"

namespace TSFCore {

// Constructors / initializers

EpetraVectorSpaceFactory::EpetraVectorSpaceFactory()
{}

EpetraVectorSpaceFactory::EpetraVectorSpaceFactory(
	const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
	)
{
	initialize(epetra_comm);
}

void EpetraVectorSpaceFactory::initialize(
	const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
	)
{
	epetra_comm_ = epetra_comm;
}

void EpetraVectorSpaceFactory::setUninitialized(
	Teuchos::RefCountPtr<const Epetra_Comm> *epetra_comm
	)
{
	if(epetra_comm) *epetra_comm = epetra_comm_;
	epetra_comm_ = Teuchos::null;
}

// Overridden from VectorSpaceFactory

Teuchos::RefCountPtr< const VectorSpace<EpetraVectorSpaceFactory::Scalar> >
EpetraVectorSpaceFactory::createVecSpc(int dim) const
{
	Teuchos::RefCountPtr<const Epetra_LocalMap>
		epetra_map = Teuchos::rcp(new Epetra_LocalMap(dim,0,*epetra_comm_));
	Teuchos::set_extra_data( epetra_comm_, &epetra_map );
	return Teuchos::rcp(new EpetraVectorSpace(epetra_map));
}

} // end namespace TSFCore
