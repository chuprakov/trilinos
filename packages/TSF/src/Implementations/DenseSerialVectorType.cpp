#include "DenseSerialVectorType.h"
#include "LAPACKGeneralMatrix.h"

#include "TSFError.h"
#include "DirectSolver.h"


using namespace TSF;

TSFVectorSpace DenseSerialVectorType::createSpace(int dimension) const
{
	return new DenseSerialVectorSpace(dimension);
}

TSFVectorSpace DenseSerialVectorType::createSpace(int dimension,
																										int nLocal,
																										int firstLocal) const
{
	return new DenseSerialVectorSpace(dimension);
}

TSFVectorSpace DenseSerialVectorType::createSpace(int dimension,
																										int nLocal,
																										const int* localIndices) const
{
	return new DenseSerialVectorSpace(dimension);
}

TSFMatrixOperator* DenseSerialVectorType::createMatrix(const TSFVectorSpace& domain,
																											 const TSFVectorSpace& range) const
{
	return new LAPACKGeneralMatrix(domain, range);
}

TSFLinearSolver DenseSerialVectorType::defaultSolver() const 
{
	return new DirectSolver();
}
