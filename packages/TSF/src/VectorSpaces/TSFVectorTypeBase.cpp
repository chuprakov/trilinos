#include "TSFVectorTypeBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFSmartPtr.h"
#include "TSFError.h"

using namespace TSF;


TSFVectorSpace TSFVectorTypeBase::createSpace(int dimension) const
{
	TSFError::raise("TSFVectorTypeBase::createSpace not "
									"implemented in base class");
	return TSFVectorSpace();
}

TSFVectorSpace TSFVectorTypeBase::createSpace(int dimension,
																							int nLocal,
																							int firstLocal) const
{
	TSFError::raise("TSFVectorTypeBase::createSpace not "
									"implemented in base class");
	return TSFVectorSpace();
}

TSFVectorSpace TSFVectorTypeBase::createSpace(int dimension,
																							int nLocal,
																							const int* localIndices) const
{
	TSFError::raise("TSFVectorTypeBase::createSpace not "
									"implemented in base class");
	return TSFVectorSpace();
}

TSFVectorSpace TSFVectorTypeBase::createSpace(int dimension,
																							int nLocal,
																							const int* localIndices,
																							int nGhost,
																							const int* ghostIndices) const
{
	TSFError::raise("TSFVectorTypeBase::createSpace not "
									"implemented in base class");
	return TSFVectorSpace();
}


TSFMatrixOperator* TSFVectorTypeBase::createMatrix(const TSFVectorSpace& domain,
																									 const TSFVectorSpace& range) const
{
	TSFError::raise("TSFVectorTypeBase::createMatrix not implemented for base class");
	return 0; // -Wall
}
