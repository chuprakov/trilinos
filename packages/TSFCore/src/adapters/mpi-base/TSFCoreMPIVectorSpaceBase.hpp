// ////////////////////////////////////////////////////////////////////////////
// TSFCoreMPIVectorSpaceBase.hpp

#ifndef TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
#define TSFCORE_MPI_VECTOR_SPACE_BASE_HPP

#include "TSFCoreMPIVectorSpaceBaseDecl.hpp"

namespace TSFCore {

///
template<class Scalar>
MPIVectorSpaceBase<Scalar>::MPIVectorSpaceBase()
	:mapCode_(-1),isInCore_(false)
{}

// Virtual methods with default implementations

template<class Scalar>
void MPIVectorSpaceBase<Scalar>::invalidateState()
{
	mapCode_  = -1;
	isInCore_ = false;
}

template<class Scalar>
Index MPIVectorSpaceBase<Scalar>::mapCode() const
{
	if(mapCode_ < 0) {
		const MPI_Comm mpiComm = this->mpiComm();
		int numProc = -1;
		MPI_Comm_size( mpiComm, &numProc );
		THROW_EXCEPTION(
			numProc!=1, std::logic_error
			,"MPIVectorSpaceBase<Scalar>::mapCode(): Error, have not implemented "
			"this method for more than one processor yet!"
			);
		mapCode_ = dim();
		isInCore_ = ( numProc == 1 );
	}
	return mapCode_;
}

// Overridden from VectorSpace

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isInCore() const
{
	return isInCore_;
}

template<class Scalar>
bool MPIVectorSpaceBase<Scalar>::isCompatible(const VectorSpace<Scalar>& vecSpc ) const
{
	if(isInCore() && vecSpc.isInCore() ) {
		return this->dim() == vecSpc.dim();
	}
	const MPIVectorSpaceBase<Scalar>
		*mpiVecSpc = dynamic_cast<const MPIVectorSpaceBase<Scalar>*>(&vecSpc);
	if(mpiVecSpc) {
		return mapCode() == mpiVecSpc->mapCode();
	}
	return false;
}

// private

template<class Scalar>
void MPIVectorSpaceBase<Scalar>::updateState() const
{
	if(mapCode_ < 0) {
		assert(numProc()==1); // ToDo: Figure out how to write a code for multi-processors
		mapCode_ = globalDim();
		isInCore_ = ( numProc() == 1 );
	}
}

} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
