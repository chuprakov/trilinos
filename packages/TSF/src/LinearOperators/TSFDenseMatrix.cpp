#include "TSFDenseMatrix.h"
#include "TSFEmptyVectorSpace.h"

using namespace TSF;

TSFDenseMatrix::TSFDenseMatrix()
	: TSFMatrixOperator(new TSFEmptyVectorSpace(), new TSFEmptyVectorSpace()), 
	isFactored_(false)
{}



