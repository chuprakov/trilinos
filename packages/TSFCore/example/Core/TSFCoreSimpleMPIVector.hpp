// /////////////////////////////////////////////////////////////////
// TSFCoreSimpleMPIVector.hpp

#ifndef TSFCORE_SIMPLE_MPI_VECTOR_HPP
#define TSFCORE_SIMPLE_MPI_VECTOR_HPP

#include "TSFCoreSimpleMPIVectorDecl.hpp"
#include "TSFExtended/src/Core/TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

template<class Scalar>
SimpleMPIVector<Scalar>::SimpleMPIVector( const MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> > &mpiSpace )
	: mpiSpace_(mpiSpace), localValues_(mpiSpace->localSubDim())
{}

// Overridden from MPIVectorBase

template<class Scalar>
MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> >
SimpleMPIVector<Scalar>::mpiSpace() const
{
	return mpiSpace_;
}

template<class Scalar>
void SimpleMPIVector<Scalar>::getLocalData( Scalar** values, ptrdiff_t* stride )
{
	*values = &localValues_[0];
	*stride = 1;
}

} // end namespace TSFCore

#endif // TSFCORE_SIMPLE_MPI_VECTOR_HPP
