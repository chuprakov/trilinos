#include "TSFMatrixOperator.h"
#include "TSFError.h"
#include "TSFVectorSpaceBase.h"


using namespace TSF;


TSFMatrixOperator::TSFMatrixOperator(const TSFVectorSpace& domain, 
																		 const TSFVectorSpace& range)
	: TSFLinearOperatorBase(domain, range), isFactored_(false)
{}

TSFMatrixOperator::~TSFMatrixOperator(){;}


void TSFMatrixOperator::setGraph(int nLocalRows, const int* bandwidth,
																 const int** columnIndices)
{
	if (requiresGraph()) 
		{
			TSFError::raise("TSFMatrixOperator::setGraph"
													 " default do-nothing version"
													 " called for a matrix that"
													 " requires a graph");
		}
}

void TSFMatrixOperator::setBandwidth(int nLocalRows, const int* bandwidth)
{
	if (requiresBandwidth()) 
		{
			TSFError::raise("TSFMatrixOperator::setBandwidth"
													 " default do-nothing version"
													 " called for a matrix that"
													 " requires bandwidth specifier");
		}	
}

void TSFMatrixOperator::setRowStructure(int globalRowIndex, int bandwidth,
																				const int* columnIndices)
{}

void TSFMatrixOperator::setRowStructure(int globalRowIndex, int bandwidth,
																				const int* columnIndices, const double* values)
{}

void TSFMatrixOperator::freezeStructure()
{}

void TSFMatrixOperator::freezeValues()
{}

void TSFMatrixOperator::setElement(int /* i */, int /* j */, const TSFReal& /* aij */)
{
	TSFError::raise("TSFMatrixOperator::setElement not avaliable "
									"for base class");
}

void TSFMatrixOperator::getILUKPreconditioner(int fillLevels,
																							int overlapFill,
																							TSFPreconditioner& rtn) const
{
	TSFError::raise("TSFMatrixOperator::getILUKPreconditioner not avaliable "
									"for base class");
}

void TSFMatrixOperator::factor()
{
	TSFError::raise("TSFMatrixOperator::factor not avaliable "
									"for base class");
}

bool TSFMatrixOperator::isFactored() const
{
	return false;
}

void TSFMatrixOperator::getRow(int row, TSFArray<int>& indices, 
                               TSFArray<TSFReal>& values) const
{
	TSFError::raise("TSFMatrixOperator::getRow not implemented "
									"for base class");
}
