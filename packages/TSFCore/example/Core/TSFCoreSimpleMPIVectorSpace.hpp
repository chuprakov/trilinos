// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSimpleMPIVectorSpace.hpp

#ifndef TSFCORE_SIMPLE_MPI_VECTOR_SPACE_HPP
#define TSFCORE_SIMPLE_MPI_VECTOR_SPACE_HPP

#include "TSFCoreSimpleMPIVectorSpaceDecl.hpp"
#include "TSFCoreSimpleMPIVector.hpp"

namespace TSFCore {

template<class Scalar>
SimpleMPIVectorSpace<Scalar>::SimpleMPIVectorSpace( MPI_Comm mpiComm, const Index localSubDim )
	:mpiComm_(mpiComm), localSubDim_(localSubDim)
{
	MPI_Comm_size( mpiComm_, &numProc_  );
	MPI_Comm_rank( mpiComm_, &procRank_ );
	globalDim_    = numProc_  * localSubDim_;
	localOffset_  = procRank_ * localSubDim_;
}

// Overridden from VectorSpece

template<class Scalar>
Index SimpleMPIVectorSpace<Scalar>::dim() const
{
	return globalDim_;
}

template<class Scalar>
MemMngPack::ref_count_ptr<Vector<Scalar> >
SimpleMPIVectorSpace<Scalar>::createMember() const
{
	return MemMngPack::rcp(new SimpleMPIVector<Scalar>(MemMngPack::rcp(this,false)));
}

template<class Scalar>
MemMngPack::ref_count_ptr< const VectorSpace<Scalar> >
SimpleMPIVectorSpace<Scalar>::clone() const
{
	return MemMngPack::rcp(new SimpleMPIVectorSpace<Scalar>(mpiComm_,localSubDim_));
}

// Overridden from MPIVectorSpaceBase

template<class Scalar>
MPI_Comm SimpleMPIVectorSpace<Scalar>::mpiComm() const
{
	return mpiComm_;
}

template<class Scalar>
Index SimpleMPIVectorSpace<Scalar>::localOffset() const
{
	return localOffset_;
}

template<class Scalar>
Index SimpleMPIVectorSpace<Scalar>::localSubDim() const
{
	return localSubDim_;
}

} // end namespace TSFCore

#endif // TSFCORE_SIMPLE_MPI_VECTOR_SPACE_HPP
