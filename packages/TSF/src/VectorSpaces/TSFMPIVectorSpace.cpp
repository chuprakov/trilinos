#include "TSFMPIVectorSpace.h"

using namespace TSF;

TSFMPIVectorSpace::TSFMPIVectorSpace()
	: TSFVectorSpaceBase()
#ifdef HAVE_MPI
	, comm_(MPI_COMM_WORLD)
#endif
{}

#ifdef HAVE_MPI
TSFMPIVectorSpace::TSFMPIVectorSpace(const MPI_Comm& comm)
	: TSFVectorSpaceBase(), comm_(comm)
{}
#endif



