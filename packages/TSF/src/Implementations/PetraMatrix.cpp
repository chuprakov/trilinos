#include "TSFDefs.h"



#define PETRA_BOOL_SUPPORTED



#include "PetraVectorSpace.h"
#include "PetraMatrix.h"
#include "TSFUtils.h"
#include "IfpackOperator.h"
#include "Ifpack_CrsRiluk.h"
#include "TSFLeftPreconditioner.h"
#include "TSFPreconditioner.h"
#include "PetraVector.h"
#include "TSFOut.h"
#include "TSFTimeMonitor.h"
#include "TSFRightPreconditioner.h"

#if HAVE_PETRA_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#endif


#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"

#include "Epetra_CrsMatrix.h"

#include "MatlabWriter.h"



using namespace TSF;


PetraMatrix::PetraMatrix(const TSFVectorSpace& domain,
												 const TSFVectorSpace& range)
	: 
	TSFMatrixOperator(domain, range), 
	rowMap_(PetraVectorSpace::getLocalMap(range)), 
	columnMap_(PetraVectorSpace::getLocalMap(domain)), 
	matrix_(0), 
	nCols_(0),
	petraComm_(0), 
	transposed_(false),
	hasCachedPreconditioner_(false),
  cachedPreconditioner_(),
  hasTranspose_(false)
{
#if HAVE_PETRA_MPI
	petraComm_ = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	petraComm_ = new Epetra_SerialComm();
#endif
}

void PetraMatrix::apply(const TSFVector& argument, 
												TSFVector& result) const 
{
	ptrCheck("apply");

	// get petra vectors from the abstract arguments
	const Epetra_Vector& in = PetraVector::getLocalValues(argument);
	Epetra_Vector& out = PetraVector::getLocalValues(result);

	// petra's mvmult is logically const but declared non-const since 
	// internal data changes. So, do a const_cast.
	Epetra_CrsMatrix* mPtr = const_cast<Epetra_CrsMatrix*>(&(*matrix_));


	// do the multiply. Set the transpose flag according to whether this
	// object represents a transpose.
	int ierr;
	{
		TSFTimeMonitor t(mvMultTimer());
		ierr = mPtr->Multiply((int) transposed_, in, out);
	}

	if (ierr != 0) 
		{
			TSFError::raise("detected error in PetraMatrix::apply, ierr=" 
											+ TSFUtils::toString(ierr));
		}

}


void PetraMatrix::applyAdjoint(const TSFVector& argument, 
															 TSFVector& result) const 
{
	ptrCheck("applyTranspose");
	// get petra vectors from the abstract arguments
	const Epetra_Vector& in = PetraVector::getLocalValues(argument);
	Epetra_Vector& out = PetraVector::getLocalValues(result);



	// petra's mvmult is logically const but declared non-const since 
	// internal data changes. So, do a const_cast.
	Epetra_CrsMatrix* mPtr = const_cast<Epetra_CrsMatrix*>(&(*matrix_));
	
	// do the multiply. Set the transpose flag according to whether this
	// object represents a transpose.

	int ierr;
	{
		TSFTimeMonitor t(mvMultTimer());
		ierr = mPtr->Multiply((int) !transposed_, in, out);
	}
	if (ierr != 0) 
		{
			TSFError::raise("detected error in PetraMatrix::applyAdjoint, ierr=" 
											+ TSFUtils::toString(ierr));
		}

}


void PetraMatrix::setBandwidth(int nLocalRows, const int* bandwidth)
{
	nCols_.resize(nLocalRows);
	for (int i=0; i<nLocalRows; i++)
		{
			nCols_[i] = *(bandwidth++);
		}
}

void PetraMatrix::setRowStructure(int globalRowIndex, int bandwidth,
																	const int* columnIndices)
{
	
			TSFArray<double> zeros(bandwidth);
			for (int i=0; i<bandwidth; i++)
				{
					zeros[i] = 0.0;
				}
            setRowStructure(globalRowIndex, bandwidth, columnIndices, &(zeros[0]));
/* 			int ierr = matrix_->InsertGlobalValues(globalRowIndex, bandwidth,  */
/* 																						 &(zeros[0]), (int*) columnIndices); */
/* 			if (ierr < 0) */
/* 				{ */
/* 					TSFError::raise("petra InsertGlobalValues row=" + TSFUtils::toString(globalRowIndex) + " failed with ierr="  */
/* 													+ TSFUtils::toString(ierr)); */
/* 				} */
	
/* 		} */
/* 	catch(exception& e) */
/* 		{ */
/* 			TSFError::trace(e, "in PetraMatrix::setRowStructure()"); */
/* 		} */
}


void PetraMatrix::setRowStructure(int globalRowIndex, int bandwidth,
                                  const int* columnIndices,
                                  const double* values)
{
	try
		{
			int ierr = matrix_->InsertGlobalValues(globalRowIndex, bandwidth, 
																						 (double*)values, (int*) columnIndices);
			if (ierr < 0)
				{
					TSFError::raise("petra InsertGlobalValues row=" + TSFUtils::toString(globalRowIndex) + " failed with ierr=" 
													+ TSFUtils::toString(ierr));
				}
	
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in PetraMatrix::setRowStructure()");
		}
}
  



void PetraMatrix::addToRow(int globalRowIndex,
													 int nCols,
													 const int* globalColumnIndices,
													 const TSFReal* a)
{
	ptrCheck("addToRow");

	// cast to non-const arrays to conform to petra
	int* indices = (int*) globalColumnIndices;
	double* values = (double*) a;

	// call petra's sumInto method
	//int ierr = matrix_->SumIntoGlobalValues(globalRowIndex, nCols, 
	//																				values, indices);
	int ierr = matrix_->SumIntoGlobalValues(globalRowIndex, nCols, 
																					values, indices);
	if (ierr < 0)
		{
			TSFError::raise("petra SumIntoGlobalValues failed with ierr=" 
											+ TSFUtils::toString(ierr));
		}
}

void PetraMatrix::setElement(int i, int j, const TSFReal& aij)
{
  ptrCheck("setElement");
  double v = aij;
  int ierr = matrix_->ReplaceGlobalValues(i, 1, &v, &j);
	if (ierr < 0)
		{
			TSFError::raise("petra setElement failed with ierr=" 
											+ TSFUtils::toString(ierr));
		}
  
}

TSFLinearOperator* PetraMatrix::getTranspose() 
{
  if (hasTranspose_) 
    {
      return &opTrp_;
    }
  hasTranspose_ = true;
  int numRowsA = range().dim();
  int numColsA = domain().dim();

  opTrp_ = new PetraMatrix(range(), domain());
  TSFSmartPtr<TSFLinearOperatorBase> Abase = opTrp_.getPtr();
  TSFMatrixOperator* ATmat = dynamic_cast<TSFMatrixOperator*>(&(*Abase));
 
  TSFArray<int> indices;
  TSFArray<double> values;
  TSFArray<TSFArray<int> > rowT(numColsA);
  TSFArray<TSFArray<double> > valT(numColsA);
  TSFArray<int> bandwidth(numColsA);
  for (int i = 0; i < numColsA; i++)
    {
      bandwidth[i] = 0;
    }
      
  for (int i = 0; i < numRowsA ; i++)
    {
      indices.resize(0);
      values.resize(0);
      getRow(i, indices, values);
      for(int j = 0; j < indices.size(); j++)
        {
          bandwidth[indices[j]]++;
          rowT[indices[j]].append(i);
          valT[indices[j]].append(values[j]);
        }
    }

  ATmat->setBandwidth(numColsA, &(bandwidth[0]));
  ATmat->freezeStructure();

  for (int i = 0; i < numColsA; i++)
    {
      ATmat->setRowStructure(i, rowT[i].size(), &(rowT[i][0]), &(valT[i][0]));
    }
/*   for (int i = 0; i < numRowsA; i++) */
/*     { */
/*       indices.resize(0); */
/*       values.resize(0); */
/*       getRow(i, indices, values); */
/*       for (int j = 0; j < indices.size(); j++) */
/*         { */
/*           ATmat->setElement(indices[j], i, values[j]); */
/*         } */
/*     } */

  ATmat->freezeValues();
  
  return &opTrp_;
}

void PetraMatrix::freezeStructure()
{
	try
		{
			int* ncols = (int*)(&(nCols_[0]));
			if (!columnMap_.isNull())
				{
					matrix_ = new Epetra_CrsMatrix(Copy, *rowMap_, ncols);
				}
			else
				{
					matrix_ = new Epetra_CrsMatrix(Copy, *rowMap_, ncols);
				}	
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in PetraMatrix::setFreezeStructure()");
		}
}

void PetraMatrix::freezeValues()
{
	int ierr = 0;
	hasCachedPreconditioner_ = false;

	if (!columnMap_.isNull())
				{
					ierr = matrix_->FillComplete(*columnMap_, *rowMap_);
				}
			else
				{
					ierr = matrix_->FillComplete();
				}
	if (ierr != 0)
		{
			TSFError::raise("petra transform failed with ierr="
											+ TSFUtils::toString(ierr));
		}
}



void PetraMatrix::zero()
{
	ptrCheck("addToRow");
	matrix_->PutScalar(0.0);
}


void PetraMatrix::print(ostream& os) const 
{
	if (matrix_.get()==0)
		{
			os << "PetraMatrix[null]";
		}
	else
		{
			os << *matrix_;
		}
}

TSFLinearOperatorBase* PetraMatrix::deepCopy() const 
{
	PetraMatrix* rtn = new PetraMatrix(domain(), range());
	rtn->petraComm_ = petraComm_;
	rtn->matrix_ = new Epetra_CrsMatrix(*matrix_);
	rtn->rowMap_ = new Epetra_Map(*rowMap_);
	rtn->columnMap_ = new Epetra_Map(*columnMap_);
	rtn->transposed_ = transposed_;

	return rtn;
}

void PetraMatrix::ptrCheck(const string& methodName) const 
{
	if (matrix_.get()==0) TSFError::raise("null matriox pointer in PetraMatrix::"
																	+ methodName);
}

void PetraMatrix::getILUKPreconditioner(int fillLevels, int overlapFill,
																				TSFPreconditioner& p) const 
{
	if (hasCachedPreconditioner_)
		{
			p = cachedPreconditioner_;
		}
	else
		{
			ptrCheck("getILUKPreconditioner");
			TSFTimeMonitor t(iluTimer());
			double relaxValue = 0.0;
			
			const Epetra_CrsGraph& matrixGraph = matrix_->Graph();
			
			Ifpack_IlukGraph* precondGraph 
			  = new Ifpack_IlukGraph(matrixGraph, fillLevels, overlapFill);
			
			petraCheck(precondGraph->ConstructFilledGraph(), 
								 "ILUK graph construct graph");
			
			Ifpack_CrsRiluk* precond 
			  = new Ifpack_CrsRiluk(*precondGraph);
			
			precond->SetRelaxValue(relaxValue);
			precond->SetRelativeThreshold(1.0);
			precond->SetAbsoluteThreshold(0.0);
			petraCheck(precond->InitValues(*matrix_), "Riluk init values");
			petraCheck(precond->Factor(), "Riluk factor");
			
			TSFLinearOperator left = new IfpackOperator(range(), domain(),
																									precond, precondGraph);
			
			p = new TSFLeftPreconditioner(left);
			cachedPreconditioner_ = p;
			hasCachedPreconditioner_ = true;
		}
}


void PetraMatrix::getILUKRightPreconditioner(int fillLevels, int overlapFill,
                                             TSFPreconditioner& p) const 
{
	if (hasCachedPreconditioner_)
		{
			p = cachedPreconditioner_;
		}
	else
		{
			ptrCheck("getILUKRightPreconditioner");
			TSFTimeMonitor t(iluTimer());
			double relaxValue = 0.0;
			
			const Epetra_CrsGraph& matrixGraph = matrix_->Graph();
			
			Ifpack_IlukGraph* precondGraph 
			  = new Ifpack_IlukGraph(matrixGraph, fillLevels, overlapFill);
			
			petraCheck(precondGraph->ConstructFilledGraph(), 
								 "ILUK graph construct graph");
			
			Ifpack_CrsRiluk* precond 
			  = new Ifpack_CrsRiluk(*precondGraph);
			
			precond->SetRelaxValue(relaxValue);
			precond->SetRelativeThreshold(1.0);
			precond->SetAbsoluteThreshold(0.0);
			petraCheck(precond->InitValues(*matrix_), "Riluk init values");
			petraCheck(precond->Factor(), "Riluk factor");
			
			TSFLinearOperator right = new IfpackOperator(range(), domain(),
                                                   precond, precondGraph);
			
			p = new TSFRightPreconditioner(right);
			cachedPreconditioner_ = p;
			hasCachedPreconditioner_ = true;
		}
}

void PetraMatrix::petraCheck(int ierr, const string& methodName) const 
{
	if (ierr < 0)
		{
			TSFError::raise("Petra error code ierr=" + TSFUtils::toString(ierr)
											+ " detected in method " + methodName);
		}
}


TSFTimer& PetraMatrix::mvMultTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("PetraMatrix mvmult");
	return *timer;
}

TSFTimer& PetraMatrix::iluTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("ILU factoring of PetraMatrix");
	return *timer;
}

Epetra_CrsMatrix* 
PetraMatrix::getConcrete(const TSFLinearOperator& A)
{
	if (A.isMatrixOperator())
		{
			const TSFSmartPtr<const TSFMatrixOperator> M = A.getMatrix();
			const PetraMatrix* pm = dynamic_cast<const PetraMatrix*>(M.get());
			if (pm==0) TSFError::raise("PetraMatrix::getConcrete bad cast");
			PetraMatrix* ptr = const_cast<PetraMatrix*>(pm);
			return ptr->matrix_.get();
		}
	else
		{
			TSFError::raise("PetraMatrix::getConcrete called for non-matrix op");
		}
	return 0; // -Wall
}

void PetraMatrix::getRow(int row, TSFArray<int>& indices, 
                                 TSFArray<TSFReal>& values) const
{
  if (isFactored_)
    {
      TSFError::raise("PetraMatrix: can't getRow since matrix is factored");
    }
  int maxRowSize = matrix_->GlobalMaxNumEntries(); 
  //cerr << "In PetraMatrix: row = " << row << endl;
  indices.resize(maxRowSize);
  values.resize(maxRowSize);
  int numEntries = 0;
  matrix_->ExtractGlobalRowCopy(row, maxRowSize, numEntries, &(values[0]), &(indices[0]));
  //cerr << "In PetraMatrix: got row " << row <<  " numEntries = "
  //  << numEntries << " maxRowSize = "<< maxRowSize << endl;
  indices.resize(numEntries);
  values.resize(numEntries);
  //cerr << "In PetraMatrix: getRow: row = " << row << " numEntries = " << numEntries << endl;
/*   for (int i = 0; i < numEntries; i++) */
/*     { */
/*       cerr << "   " << indices[i] << "  " << values[i] << endl; */
/*     } */
}











