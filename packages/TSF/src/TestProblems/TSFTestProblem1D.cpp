#include "TSFTestProblem1D.h"

#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFMPI.h"
#include <cmath>

using namespace TSF;


TSFTestProblem1D::TSFTestProblem1D(TSFReal a, TSFReal b, int nGlobal, 
																	 TSFMatrixOperator* matrix)
	: TSFDefaultMatrixProblem(matrix),
		h_((b-a)/((TSFReal) (nGlobal-1))),
		a_(a),
		b_(b),
		nGlobal_(nGlobal),
		nLocal_(0),
		myLowestRow_(0)
{
	int nProc = TSFMPI::getNProc();
	int myRank = TSFMPI::getRank();

	int rowsPerProc = (int) floor(((double) nGlobal)/((double) nProc));
	myLowestRow_ = myRank*rowsPerProc;
	int myHighestRow = rowsPerProc + myLowestRow_;
	if (myRank==nProc-1)
		{
			myHighestRow = nGlobal_;
		}
	nLocal_ = myHighestRow - myLowestRow_;
}

TSFReal TSFTestProblem1D::getX(int globalIndex) const 
{
	return a_ + (b_-a_)*((TSFReal) globalIndex)/((TSFReal)(nGlobal_-1));
}


void TSFTestProblem1D::getRowValues(int row, TSFArray<int>& indices,
																		TSFArray<TSFReal>& values) const 
{
	indices.resize(3);
	values.resize(3);

	if (row==0)
		{
			indices[0] = row;
			indices[1] = row+1;
			indices[2] = row+2;
			values[0] = leftBCZerothOrderTerm() - 1.5*leftBCFirstOrderTerm()/h_;
			values[1] = 2.0*leftBCFirstOrderTerm()/h_;
			values[2] = -0.5*leftBCFirstOrderTerm()/h_;
		}
	else if (row==(nGlobal_-1))
		{
			indices[0] = row;
			indices[1] = row-1;
			indices[2] = row-2;
			values[0] = -rightBCZerothOrderTerm() + 1.5*rightBCFirstOrderTerm()/h_;
			values[1] = -2.0*rightBCFirstOrderTerm()/h_;
			values[2] = 0.5*rightBCFirstOrderTerm()/h_;
		}
	else
		{
			indices[0] = row-1;
			indices[1] = row;
			indices[2] = row+1;

			TSFReal x0 = getX(row);
			TSFReal xPlus = x0 + h_;
			TSFReal xMinus = x0 - h_;
			values[0] = (secondOrderTerm(xMinus)/h_ - firstOrderTerm(xMinus)/2.0)/h_;
			values[1] = zerothOrderTerm(x0) - 2.0*secondOrderTerm(x0)/h_/h_;
			values[2] = (secondOrderTerm(xPlus)/h_ + firstOrderTerm(xPlus)/2.0)/h_;
		}

}


void TSFTestProblem1D::getRowStruct(int row, TSFNonDupTSFArray<int>& indices) const 
{
	indices.resize(3);

	if (row==0)
		{
			indices[0] = row;
			indices[1] = row+1;
			indices[2] = row+2;
		}
	else if (row==(nGlobal_-1))
		{
			indices[0] = row;
			indices[1] = row-1;
			indices[2] = row-2;
		}
	else
		{
			indices[0] = row-1;
			indices[1] = row;
			indices[2] = row+1;
		}
}

TSFReal TSFTestProblem1D::getRHSValue(int row) const 
{
	if (row==0)
		{
			return leftBCForcingTerm();
		}
	else if (row==(nGlobal_-1))
		{
			return rightBCForcingTerm();
		}
	else
		{
			return forcingTerm(getX(row));
		}
}

TSFReal TSFTestProblem1D::getSolutionValue(int row) const 
{
	return solution(getX(row));
}


		
