#include "TSFMPIVector.h"
#include "TSFError.h"
#include "TSFMPIVectorSpace.h"


using namespace TSF;

TSFMPIVector::TSFMPIVector(const TSFVectorSpace& space)
	: TSFInCoreVector(space)
#if HAVE_MPI
	, comm_()
#endif
{
#if HAVE_MPI
	const TSFMPIVectorSpace* mpiSpace 
		= dynamic_cast<const TSFMPIVectorSpace*>(space.ptr());
	if (mpiSpace==0) 
		{
			TSFError::raise("TSFMPIVector attempted to cast non-mpi vector space "
											"to TSFMPIVectorSpace");
		}
	comm_ = mpiSpace->getMPIComm();
#endif
}

TSFReal TSFMPIVector::reduce(const TSFReal& localValue, TSFReductionOp op) const
{
#if HAVE_MPI
	TSFTimeMonitor timer(commTimer());

	MPI_Op mpiop =  MPI_SUM;
	switch(op)
		{
		case TSFSumReduce: 
			mpiop = MPI_SUM;
			break;
		case TSFMaxReduce: 
			mpiop = MPI_MAX;
			break;
		}
	double reducedValue;
	int ierr = MPI_Allreduce((void*) &localValue, (void*) &reducedValue, 1, 
													 MPI_DOUBLE, mpiop, comm_);
	if (ierr < 0) TSFError::raise("error in TSFMPIVector::reduce");

	return reducedValue;
												
#else
	return localValue;
#endif
}

TSFTimer& TSFMPIVector::commTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("TSFMPIVector comm");
	return *timer;
}
