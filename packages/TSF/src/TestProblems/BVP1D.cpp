#include "BVP1D.h"

#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include <math.h>

using namespace TSF;

BVP1D::BVP1D(int n, TSFReal a, TSFReal b,
						 TSFMatrixOperator* matrix)
	: TSFDefaultMatrixProblem(matrix), n_(n), a_(a), b_(b),
		h_((b-a)/((TSFReal) n-1))
{
	
	int nProc = TSFMPI::getNProc();
	int myPid = TSFMPI::getRank();

	int nRows =	n_;
	int rowsPerProc = (int) floor(((double) nRows)/((double) nProc));
	int myLowestRow = myPid*rowsPerProc;
	int myHighestRow = rowsPerProc + myLowestRow;
	if (myPid==nProc-1)
		{
			myHighestRow = nRows;
		}

	nLocal_ = myHighestRow - myLowestRow;

	lowestLocalRow_ = myLowestRow;

	TSFOut::printf("low=%d high=%d\n", lowestLocalRow_, myHighestRow);
}

TSFSmartPtr<TSFArray<int> > BVP1D::formUpdateList() const
{
	TSFSmartPtr<TSFArray<int> > rtn = new TSFArray<int>(nLocal_);
	for (int i=0; i<nLocal_; i++)
		{
			(*rtn)[i] = lowestLocalRow_ + i;
		}
	return rtn;
}

int BVP1D::getRowBandwidth(int row) const 
{
	if (row==0 || row==n_-1) return 1;
	return 3;
}

void BVP1D::getRowValues(int row, TSFArray<int>& indices, 
												 TSFArray<TSFReal>& values) const
{
	indices.resize(getRowBandwidth(row));
	values.resize(getRowBandwidth(row));

	TSFReal w = 1.0/h()/h();



	if (row==0 || row==n_-1) 
		{
			indices[0] = row;
			values[0] = w;
		}
	else
		{
			TSFReal xMinus = (getX(row-1) + getX(row))/2.0;
			TSFReal xPlus = (getX(row+1) + getX(row))/2.0;
			indices[0] = row-1;
			indices[1] = row;
			indices[2] = row+1;
			values[0] = w*alpha(xMinus);
			values[1] = -w*(alpha(xMinus) + alpha(xPlus)) + beta(getX(row));
			values[2] = w*alpha(xPlus);
		}
}

TSFReal BVP1D::getRHSValue(int row) const 
{
	if (row==0)
		{
			return uLeft()/h()/h();
		}
	else if (row==n_-1)
		{
			return uRight()/h()/h();
		}
	else 
		{
			return gamma(getX(row));
		}
}

TSFReal BVP1D::getSolutionValue(int row) const 
{
	return solution(getX(row));
}



