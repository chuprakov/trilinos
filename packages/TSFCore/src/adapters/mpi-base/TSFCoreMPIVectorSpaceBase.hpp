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
		const Index localSubDim = this->localSubDim(); 
		int numProc = 1;
		int procRank = 0;
#ifdef RTOp_USE_MPI
		if( mpiComm != MPI_COMM_NULL ) {
			MPI_Comm_size( mpiComm, &numProc );
			MPI_Comm_rank( mpiComm, &procRank );
		}
		if(numProc > 1) {
			//
			// Here we will make a map code out of just the local
			// sub-dimension on each processor.  If each processor
			// has the same number of local elements, then the maps
			// will be the same and this is all you need for
			// RTOp compatibility unless the operations are not
			// coordinate invariant.  I will work on this issue
			// if it becomes a problem.
			//
			Index localCode = localSubDim % procRank + localSubDim;
			int *dummy = &localCode; // For now make sure that Index is int so we can use MPI_INT
			MPI_Allreduce(
				&localCode            // sendbuf
				,&mapCode_            // recvbuf
				,1                    // count
				,MPI_INT              // datatype (ToDo: use traits class on Index to make more general)
				,MPI_SUM              // op
				,mpiComm              // comm
				);
			isInCore_ = false;
		}
		else {
#endif	
			mapCode_ = localSubDim;
			isInCore_ = true;
#ifdef RTOp_USE_MPI
		}
#endif
	}
    else {
		mapCode_  = -1;     // Uninitialized!
		isInCore_ = false;
	}
}
	
} // end namespace TSFCoreo

#endif // TSFCORE_MPI_VECTOR_SPACE_BASE_HPP
