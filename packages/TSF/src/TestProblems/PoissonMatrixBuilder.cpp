#include "PoissonMatrixBuilder.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"



#include <math.h>

using namespace TSF;

PoissonMatrixBuilder::PoissonMatrixBuilder(int nGlobal, const PComm& comm)
	: nGlobal_(nGlobal), 
	updates_(0), 
	columns_(0), 
	comm_(comm), 
	h_(0.0), 
	nProc_(1),
	myPid_(0), 
	myNRows_(0)
 
{
	nProc_ = comm_.getNProc();
	myPid_ = comm_.getRank();
	h_ = 1.0/((double) (nGlobal_-1));

	/* find the range of rows for which this processor is responsible */
	int nRows = nGlobal_;
	int rowsPerProc = (int) floor(((double) nRows)/((double) nProc_));
	int myLowestRow = myPid_*rowsPerProc;
	int myHighestRow = rowsPerProc + myLowestRow;
	if (myPid_==nProc_-1)
		{
			myHighestRow = nRows;
		}
	myNRows_ = myHighestRow - myLowestRow;
	
	PMachine::printf("proc %d is responsible for rows %d-%d of %d",
									 myPid_, myLowestRow, myHighestRow-1, nRows);

	/* fill in the update list and graph */
	updates_ = new TSFArray<int>(myNRows_);
	columns_.resize(myNRows_);

	for (int i=0; i<myNRows_; i++)
		{
			int row = myLowestRow + i;
			(*updates_)[i] = row;
			getRowStruct(row, columns_[i]);
		}
}

void PoissonMatrixBuilder::buildSystem(TSFLinearOperator& op, 
																					 TSFVector& rhs,
																					 TSFVector& ans) const 
{
	TSFLinearOperatorBase* ptr = &(*op.getPtr());
	MatrixOperatorBase* A = dynamic_cast<MatrixOperatorBase*>(ptr);
	
	buildMatrix(A);

	rhs = op.range().createMember();
	ans = op.domain().createMember();

	buildRHS(rhs);
	buildAnswer(ans);
}

void PoissonMatrixBuilder::buildMatrix(MatrixOperatorBase* A) const 
{
	/* set the matrix communicator to our communicator */
	A->setCommunicator(comm_);

	/* set the row and column spaces of the matrix */
	A->setRowMap(nGlobal_, updates_);
	A->setColumnMap(nGlobal_, updates_);
	
	/* if the graph is needed for this matrix type, set it */
	if (A->requiresGraph())
		{
			A->setGraph(columns_);
		}

	/* tell the matrix that we're done with configuration */
	A->freezeStructure();


	/* set all elements to zero */
	A->zero();


	/* load the matrix */
	for (int i=0; i<updates_->length(); i++)
		{
			int row = (*updates_)[i];
			TSFArray<int> indices;
			TSFArray<TSFReal> values;
			getRow(row, indices, values);
			A->addToRow(row, indices, values);
		}

	A->freezeValues();
}


void PoissonMatrixBuilder::buildRHS(TSFVector& rhs) const 
{
	for (int i=0; i<updates_->length(); i++)
		{
			int row = (*updates_)[i];
			PMachine::printf("setting element %d to %g", row, getRHS(row));
			rhs.setElement(row, getRHS(row));
		}
}

void PoissonMatrixBuilder::buildAnswer(TSFVector& ans) const 
{
	for (int i=0; i<updates_->length(); i++)
		{
			int row = (*updates_)[i];
			ans.setElement(row, getSoln(row));
		}
}



double PoissonMatrixBuilder::getRHS(int row) const 
{
	double x = row*h_;
	
	if (row==0 || row==nGlobal_-1) return 0.0;
	return 2.0-6.0*x;
}

double PoissonMatrixBuilder::getSoln(int row) const 
{
	double x = row*h_;
	return x*x*(1.0-x);
}

void PoissonMatrixBuilder::getRowStruct(int row, 
																						NonDupTSFArray<int>& indices) const 
{
	if (row==0 || row==(nGlobal_-1))
		{
			indices.append(Int(row));
		}
	else
		{
			indices.append(Int(row-1));
			indices.append(Int(row));
			indices.append(Int(row+1));
		}
}

void PoissonMatrixBuilder::getRow(int row, 
																			TSFArray<int>& indices, 
																			TSFArray<TSFReal>& values) const 
{
	if (row==0 || row==(nGlobal_-1))
		{
			indices.append(Int(row));
			values.append(1.0/h_/h_);
		}
	else
		{
			indices.append(Int(row-1));
			values.append(1.0/h_/h_);
			indices.append(Int(row));
			values.append(-2.0/h_/h_);
			indices.append(Int(row+1));
			values.append(1.0/h_/h_);
		}
}
