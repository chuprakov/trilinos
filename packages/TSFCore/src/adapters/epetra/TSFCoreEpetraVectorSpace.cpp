// ////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpace.cpp

// Define this to use the optimized EpetraMultiVector implementation!
#define TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR

#include <assert.h>

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVectorSpaceFactory.hpp"
#include "TSFCoreEpetraVector.hpp"
#ifdef TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR
#  include "TSFCoreEpetraMultiVector.hpp"
#endif
#include "Teuchos_TestForException.hpp"
#include "dynamic_cast_verbose.hpp"

#include "Epetra_Comm.h"
#ifdef RTOp_USE_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Map.h"

namespace TSFCore {

EpetraVectorSpace::EpetraVectorSpace()
{
	updateState();
}

EpetraVectorSpace::EpetraVectorSpace(
	const Teuchos::RefCountPtr<const Epetra_Map>  &epetra_map
	)
{
	initialize(epetra_map);
}

void EpetraVectorSpace::initialize(
	const Teuchos::RefCountPtr<const Epetra_Map>  &epetra_map
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_map.get(), std::invalid_argument, "EpetraVectorSpace::initialize(...): Error!" );
#endif // _DEBUG
	epetra_map_  = epetra_map;
#ifdef RTOp_USE_MPI
	const Epetra_MpiComm
		*epetra_mpi_comm = dynamic_cast<const Epetra_MpiComm*>(&epetra_map_->Comm());
	if(epetra_mpi_comm) {
//		std::cout << "EpetraVectorSpace::initialize(...): Using an Epetra_MpiComm!\n";
		mpiComm_ = epetra_mpi_comm->Comm();
#ifdef _DEBUG
		TEST_FOR_EXCEPTION(
			mpiComm_ == MPI_COMM_NULL, std::logic_error
			,"EpetraVectorSpace::initialize(...), Error, if using Epetra_MpiComm then "
			"the associated MPI_Comm object can not be MPI_COMM_NULL!"
			);
#endif // _DEBUG
//		std::cout << "EpetraVectorSpace::initialize(...): mpiComm_ = " << mpiComm_ << std::endl;
	}
	else {
		std::cout << "EpetraVectorSpace::initialize(...): Not using an Epetra_MpiComm!\n";
		mpiComm_ = MPI_COMM_NULL;
	}
#else // RTOp_USE_MPI
//	std::cout << "EpetraVectorSpace::initialize(...): Not using an Epetra_MpiComm!\n";
	mpiComm_ = MPI_COMM_NULL;
#endif // RTOp_USE_MPI
	localSubDim_ = epetra_map->NumMyElements();
	smallVecSpcFcty_ = Teuchos::rcp(
		new EpetraVectorSpaceFactory(
			Teuchos::rcp(&epetra_map->Comm(),false)
			)
		);
	updateState();
}

void EpetraVectorSpace::setUninitialized(
	Teuchos::RefCountPtr<const Epetra_Map> *epetra_map
	)
{
	if(epetra_map) *epetra_map = epetra_map_;
	epetra_map_ = Teuchos::null;
	mpiComm_ = MPI_COMM_NULL;
	localSubDim_ = 0; // Flag that this is uninitialized
	smallVecSpcFcty_ = Teuchos::null;
	updateState();
}

// Overridden from VectorSpace

Index EpetraVectorSpace::dim() const
{
	return epetra_map_.get() ? epetra_map_->NumGlobalElements() : 0;
}

Teuchos::RefCountPtr<Vector<EpetraVectorSpace::Scalar> >
EpetraVectorSpace::createMember() const
{
	return Teuchos::rcp(
    new EpetraVector(
      Teuchos::rcp(new Epetra_Vector(*epetra_map_,false))
      ,Teuchos::rcp(this,false)
      )
    );
}

Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> >
EpetraVectorSpace::smallVecSpcFcty() const
{
	return smallVecSpcFcty_;
}

Teuchos::RefCountPtr<MultiVector<EpetraVectorSpace::Scalar> > 
EpetraVectorSpace::createMembers(int numMembers) const
{
#ifdef TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR
	// Use specialized Epetra_MultiVector implementation
	return Teuchos::rcp(
		new EpetraMultiVector(
			Teuchos::rcp(new Epetra_MultiVector(*epetra_map_,numMembers,false))
			,Teuchos::rcp(this,false)
			)
		);
#else
	// Use default MultiVector implementation
	return VectorSpace<Scalar>::createMembers(numMembers);
#endif
}

Teuchos::RefCountPtr< const VectorSpace<EpetraVectorSpace::Scalar> >
EpetraVectorSpace::clone() const
{
	return Teuchos::rcp( new EpetraVectorSpace(epetra_map_) );
}

// Overriddend from MPIVectorSpaceBase

MPI_Comm EpetraVectorSpace::mpiComm() const
{
	return mpiComm_;
}

Index EpetraVectorSpace::localSubDim() const
{
	return localSubDim_;
}

} // end namespace TSFCore
