#include "TSFDefaultMatrixProblem.h"

#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"


using namespace TSF;



TSFDefaultMatrixProblem::TSFDefaultMatrixProblem(TSFMatrixOperator* matrix)
	: TSFMatrixProblem(matrix)
{}


TSFSmartPtr<TSFArray<int> > TSFDefaultMatrixProblem::formBandwidth() const 
{
	TSFSmartPtr<TSFArray<int> > rtn = new TSFArray<int>(nLocalRows());

	TSFArray<int>& b = *rtn;
	
	for (int i=0; i<nLocalRows(); i++)
		{
			int row = lowestLocalRow() + i;
			b[i] = getRowBandwidth(row);
		}

	return rtn;
}
void TSFDefaultMatrixProblem::fillMatrix() const 
{
	for (int i=0; i<nLocalRows(); i++)
		{
			int row = lowestLocalRow() + i;
			TSFArray<int> indices;
			TSFArray<TSFReal> values;
			getRowValues(row, indices, values);
			matrix_->setRowStructure(row, indices.size(), &(indices[0]));
			matrix_->addToRow(row, indices.size(), &(indices[0]), &(values[0]));
		}
}

void TSFDefaultMatrixProblem::fillRHS(TSFVector& rhs) const 
{
	for (int i=0; i<nLocalRows(); i++)
		{
			int row = lowestLocalRow() + i;
			rhs[row] = getRHSValue(row);
		}	
}

void TSFDefaultMatrixProblem::fillKnownSolution(TSFVector& soln) const 
{
	for (int i=0; i<nLocalRows(); i++)
		{
			int row = lowestLocalRow() + i;
			soln[row] = getSolutionValue(row);
		}	
}

