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
Index MPIVectorSpaceBase<Scalar>::mapCode() const
{
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

// protected

template<class Scalar>
void MPIVectorSpaceBase<Scalar>::updateState()
{
	if( this->localSubDim() > 0 ) {
		const MPI_Comm mpiComm = this->mpiComm();
		int numProc = 1;
#ifdef RTOp_USE_MPI
		MPI_Comm_size( mpiComm, &numProc );
#endif	
		TEST_FOR_EXCEPTION(
			numProc!=1, std::logic_error
			,"MPIVectorSpaceBase<Scalar>::updateState(): Error, have not implemented "
			"this method for more than one processor yet!"
			);
		mapCode_ = this->dim();
		isInCore_ = ( numProc == 1 );
	}
    else {
		mapCode_  = -1;     // Uninitialized!
		isInCore_ = false;
	}
}
	
} // end namespace TSFCore

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
