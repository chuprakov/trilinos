#include "LAPACKMatrixBase.h"


using namespace TSF;

LAPACKMatrixBase::LAPACKMatrixBase()
	: nRows_(0), nCols_(0), data_()
{}

void LAPACKMatrixBase::configureDense(int nRows, int nCols)
{
	nRows_ = nRows;
	nCols_ = nCols;
	domain_ = new DenseSerialVectorSpace(nCols);
	range_ = new DenseSerialVectorSpace(nRows);
	
}

