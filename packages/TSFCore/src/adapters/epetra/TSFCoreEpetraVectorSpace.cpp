// ////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpace.cpp

#include <assert.h>

#include "TSFCoreEpetraVectorSpace.hpp"
#include "TSFCoreEpetraVector.hpp"
//#include "TSFCoreEpetraMultiVector.hpp"
#include "ThrowException.hpp"
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
	const MemMngPack::ref_count_ptr<const Epetra_BlockMap>  &epetra_map
	)
	:mpiComm_(MPI_COMM_NULL),localOffset_(-1),localSubDim_(-1)
{
	initialize(epetra_map);
}

void EpetraVectorSpace::initialize(
	const MemMngPack::ref_count_ptr<const Epetra_BlockMap>  &epetra_map
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	THROW_EXCEPTION( !epetra_map.get(), std::invalid_argument, "EpetraVectorSpace::initialize(...): Error!" );
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

MemMngPack::ref_count_ptr<const Epetra_BlockMap>
EpetraVectorSpace::setUninitialized()
{
	MemMngPack::ref_count_ptr<const Epetra_BlockMap> tmp = epetra_map_;
	epetra_map_ = MemMngPack::null;
	return tmp;
}

// Overridden from VectorSpace

Index EpetraVectorSpace::dim() const
{
	return epetra_map_.get() ? epetra_map_->NumGlobalElements() : 0;
}

MemMngPack::ref_count_ptr<Vector<EpetraVectorSpace::Scalar> >
EpetraVectorSpace::createMember() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp(new EpetraVector(mmp::rcp(new Epetra_Vector(*epetra_map_))));
}

MemMngPack::ref_count_ptr<MultiVector<EpetraVectorSpace::Scalar> > 
EpetraVectorSpace::createMembers(int numMembers) const
{
	namespace mmp = MemMngPack;
	return VectorSpace<Scalar>::createMembers(numMembers); // Todo: use the below code!
/*
	return mmp::rcp(
		new EpetraMultiVector(
			mmp::rcp(new Epetra_MultiVector(*epetra_map_,numMembers))  // epetra_multi_vec
			,mmp::rcp(this,false)                                      // epetra_vec_spc
			)
		);
*/
}

MemMngPack::ref_count_ptr< const VectorSpace<EpetraVectorSpace::Scalar> >
EpetraVectorSpace::clone() const
{
	namespace mmp = MemMngPack;
	return mmp::rcp( new EpetraVectorSpace(epetra_map_) );
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
