#include "TSFOneDTestProblem.h"

#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"

TSFOneDTestProblem::TSFOneDTestProblem(TSFReal a, TSFReal b, int nGlobal, 
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
	if (myPid==nProc-1)
		{
			myHighestRow = nRows;
		}
	nLocal_ = myHighestRow - myLowestRow_;
}
		
