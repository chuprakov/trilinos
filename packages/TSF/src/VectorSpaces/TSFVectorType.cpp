#include "TSFVectorType.h"
#include "TSFVectorSpaceBase.h"

#include <typeinfo>

using namespace TSF;

TSFVectorType::TSFVectorType()
	: ptr_(0)
{}

TSFVectorType::TSFVectorType(TSFVectorTypeBase* ptr)
	: ptr_(ptr)
{}

TSFVectorSpace TSFVectorType::createSpace(int dimension) const
{
	return ptr_->createSpace(dimension);
}

TSFVectorSpace TSFVectorType::createSpace(int dimension,
																					int nLocal,
																					int firstLocal) const
{
	return ptr_->createSpace(dimension, nLocal, firstLocal);
}

TSFVectorSpace TSFVectorType::createSpace(int dimension,
																					int nLocal,
																					const int* localIndices) const
{
	return ptr_->createSpace(dimension, nLocal, localIndices);
}

TSFVectorSpace TSFVectorType::createSpace(int dimension,
																					int nLocal,
																					const int* localIndices,
																					int nGhost,
																					const int* ghostIndices) const
{
	if (nGhost==0) return createSpace(dimension, nLocal, localIndices);

	return ptr_->createSpace(dimension, nLocal, localIndices,
													 nGhost, ghostIndices);
}

TSFMatrixOperator* TSFVectorType::createMatrix(const TSFVectorSpace& domain,
																							 const TSFVectorSpace& range) const
{
	return ptr_->createMatrix(domain, range);
}

bool TSFVectorType::operator==(const TSFVectorType& other) const
{
	return (typeid(*other.ptr_) == typeid(*ptr_));
}

ostream& TSFVectorType::print(ostream& os) const 
{
	return ptr_->print(os);
}

TSFLinearSolver TSFVectorType::defaultSolver() const 
{
	return ptr_->defaultSolver();
}





