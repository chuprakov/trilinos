#include "LAPACKGeneralMatrix.h"

#include "TSFBlas.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFUtils.h"
#include "DenseSerialVectorSpace.h"
#include "TSFSerialVector.h"

using namespace TSF;


LAPACKGeneralMatrix::LAPACKGeneralMatrix(const TSFVectorSpace& domain,
																				 const TSFVectorSpace& range)
	: TSFMatrixOperator(domain, range), 
	  data_(range.dim()*domain.dim()),
	  iPiv_(),
	  is_factored_(false),
	  nRows_(range.dim()), 
	  nCols_(domain.dim())
{}






void LAPACKGeneralMatrix::apply(const TSFVector& in,
																TSFVector& out) const
{
	mvMult(false, in, out);
}

void LAPACKGeneralMatrix::applyAdjoint(const TSFVector& in,
																			 TSFVector& out) const
{
	mvMult(true, in, out);
}

void LAPACKGeneralMatrix::applyInverse(const TSFVector& in,
																			 TSFVector& out) const
{
	solve(false, in, out);
}

void LAPACKGeneralMatrix::applyInverseAdjoint(const TSFVector& in,
																							TSFVector& out) const
{
	solve(true, in, out);
}



void LAPACKGeneralMatrix::mvMult(bool transpose, const TSFVector& in,
																 TSFVector& out) const
{
	int one = 1;
	TSFReal onePointZero = 1.0;
	TSFReal zero = 0.0;

	const DenseSerialVector& vIn = TSFSerialVector::getConcrete(in);
	DenseSerialVector& vOut = TSFSerialVector::getConcrete(out);
    //cout << " vIn befor\n";
    //cout << vIn;
    //vIn.print(cout);
    //cout << endl;
	const TSFReal* inPtr = &(vIn[0]);
	TSFReal* outPtr = &(vOut[0]);
	
	/* set the LAPACK transpose flag = "N" for no transpose, "T" for transpose */
	char transFlag='N';
	if (transpose) transFlag='T';
	
	// RAB & ADP : 7/10/2002 : We have fixed this!
	TSFBlas<TSFReal>::gemv(&transFlag, &nRows_, &nCols_, &onePointZero, 
						   &(data_[0]),
						   &nRows_, inPtr, &one, &zero, 
						   outPtr, &one);
}

void LAPACKGeneralMatrix::solve(bool transpose, const TSFVector& in,
																TSFVector& out) const
{
	int one = 1;
	int info = 0;

	const DenseSerialVector& vIn = TSFSerialVector::getConcrete(in);
	DenseSerialVector& vOut = TSFSerialVector::getConcrete(out);

	const TSFReal* inPtr = &(vIn[0]);
	TSFReal* outPtr = &(vOut[0]);

	/* LAPACK overwrites the input vector argument. We copy the input 
	 * vector into the output vector, and then pass the output vector
	 * to the backsolve routine. */
	TSFBlas<TSFReal>::copy(&nRows_, inPtr, &one, outPtr, &one);

	/* factor if we haven't already done so */
	if (!is_factored_)
		{
			const_cast<LAPACKGeneralMatrix*>(this)->factor();
		}
	TSFReal* dataPtr = const_cast<TSFReal*>(&(factor_data_[0]));

	/* set the LAPACK transpose flag = "N" for no transpose, "T" for transpose */
	char transFlag='N';
	if (transpose) transFlag='T';

	/* backsolve */
	int* pivPtr = const_cast<int*>(&(iPiv_[0]));
	TSFBlas<TSFReal>::getrs(&transFlag, &nRows_, &one, dataPtr,
													&nRows_, pivPtr, outPtr, &nRows_, &info);

	if (info != 0)
		{
			TSFError::raise("LAPACKGeneralMatrix backsolve failed with error code"
											+ TSFUtils::toString(info));
		}
}

void LAPACKGeneralMatrix::addToRow(int globalRowIndex,
																	 int nCols,
																	 const int* globalColumnIndices,
																	 const TSFReal* a)
{
	is_factored_ = false;
	for (int i=0; i<nCols; i++)
		{
			data_[nRows_*globalColumnIndices[i] + globalRowIndex] += a[i];
		}
}

void LAPACKGeneralMatrix::setElement(int i, int j, const TSFReal& aij)
{
	is_factored_ = false;
	data_[nRows_*j + i] = aij;
}

void LAPACKGeneralMatrix::zero()
{
	is_factored_ = false;
	data_.zero();
}

void LAPACKGeneralMatrix::factor() 
{
	int info = 0;

	iPiv_.resize(nRows_);

	int* pivPtr = const_cast<int*>(&(iPiv_[0]));

//	factor_data_.resize(data_.length());
	factor_data_ = data_;

	TSFReal* dataPtr = const_cast<TSFReal*>(&(factor_data_[0]));

	TSFBlas<TSFReal>::getrf(&nRows_, &nCols_, dataPtr,
													&nRows_, pivPtr, &info);

	if (info != 0)
		{
			TSFError::raise("LAPACKGeneralMatrix factor failed with error code"
											+ TSFUtils::toString(info));
		}
	is_factored_ = true;
}

void LAPACKGeneralMatrix::print(ostream& os) const
{
	os << "LAPACK " << nRows_ << "-by-" << nCols_ << " matrix: " << endl;
	os << "[";
	for (int i=0; i<nRows_; i++)
		{
			os << "[";
			for (int j=0; j<nCols_; j++)
				{
					os << data_[i + nRows_*j];
					if (j < nCols_-1) os << ", ";
				}
			os << "]";
			if (i < nRows_-1) os << ", ";
		}
	os << "]";
}
