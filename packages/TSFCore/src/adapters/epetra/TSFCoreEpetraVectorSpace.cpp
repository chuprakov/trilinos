// ////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpace.cpp

#include <assert.h>

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
//#include "TSFCoreEpetraMultiVector.hpp"
#include "Teuchos_TestForException.hpp"
#include "dynamic_cast_verbose.hpp"

#include "Epetra_Comm.h"
#ifdef PETRA_COMM_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"


namespace TSFCore {

EpetraVectorSpace::EpetraVectorSpace()
	:mpiComm_(MPI_COMM_NULL),localOffset_(-1),localSubDim_(-1)
{}

EpetraVectorSpace::EpetraVectorSpace(
	const Teuchos::RefCountPtr<const Epetra_BlockMap>  &epetra_map
	)
	:mpiComm_(MPI_COMM_NULL),localOffset_(-1),localSubDim_(-1)
{
	initialize(epetra_map);
}

void EpetraVectorSpace::initialize(
	const Teuchos::RefCountPtr<const Epetra_BlockMap>  &epetra_map
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !epetra_map.get(), std::invalid_argument, "EpetraVectorSpace::initialize(...): Error!" );
#endif
	epetra_map_  = epetra_map;
#ifdef PETRA_COMM_MPI
	mpiComm_ = dyn_cast<const Epetra_MpiComm>(epetra_map_->Comm()).Comm();
#else
	mpiComm_ = MPI_COMM_NULL;
#endif
	localOffset_ = epetra_map->MinMyGID() - epetra_map->IndexBase();
	localSubDim_ = epetra_map->NumMyElements();
	MPIVectorSpaceBase<Scalar>::invalidateState();
}

Teuchos::RefCountPtr<const Epetra_BlockMap>
EpetraVectorSpace::setUninitialized()
{
	Teuchos::RefCountPtr<const Epetra_BlockMap> tmp = epetra_map_;
	epetra_map_ = Teuchos::null;
	return tmp;
}

// Overridden from VectorSpace

Index EpetraVectorSpace::dim() const
{
	return epetra_map_.get() ? epetra_map_->NumGlobalElements() : 0;
}

Teuchos::RefCountPtr<Vector<EpetraVectorSpace::Scalar> >
EpetraVectorSpace::createMember() const
{
	return Teuchos::rcp(new EpetraVector(Teuchos::rcp(new Epetra_Vector(*epetra_map_))));
}

Teuchos::RefCountPtr<MultiVector<EpetraVectorSpace::Scalar> > 
EpetraVectorSpace::createMembers(int numMembers) const
{
	return VectorSpace<Scalar>::createMembers(numMembers); // Todo: use the below code!
/*
	return Teuchos::rcp(
		new EpetraMultiVector(
			Teuchos::rcp(new Epetra_MultiVector(*epetra_map_,numMembers))  // epetra_multi_vec
			,Teuchos::rcp(this,false)                                      // epetra_vec_spc
			)
		);
*/
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

Index EpetraVectorSpace::localOffset() const
{
	return localOffset_;
}

Index EpetraVectorSpace::localSubDim() const
{
	return localSubDim_;
}

} // end namespace TSFCore
