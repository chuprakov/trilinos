#include "TSFMultiVectorOperator.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFUtils.h"
#include "DenseSerialVectorSpace.h"
#include "TSFSerialVector.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;

TSFMultiVectorOperator::TSFMultiVectorOperator()
	: TSFLinearOperatorBase(), isVertical_(false), vectors_(0)
{}

TSFMultiVectorOperator::TSFMultiVectorOperator(const TSFVectorSpace& space,
																							 int numVectors,
																							 bool isVertical)
	: TSFLinearOperatorBase(), isVertical_(isVertical), vectors_(numVectors)
{
	if (isVertical_)
		{
			range_ = space;
			domain_ = new DenseSerialVectorSpace(numVectors);
		}
	else
		{
			domain_ = space;
			range_ = new DenseSerialVectorSpace(numVectors);
		}

	for (int i=0; i<numVectors; i++)
		{
			vectors_[i] = space.createMember();
		}
}

void TSFMultiVectorOperator::apply(const TSFVector& in, 
																	 TSFVector& out) const
{
	if (isVertical_)
		{
			verticalMVMult(in, out);
		}
	else
		{
			horizontalMVMult(in, out);
		}
}

void TSFMultiVectorOperator::applyAdjoint(const TSFVector& in, 
																					TSFVector& out) const
{
	if (isVertical_)
		{
			horizontalMVMult(in, out);
		}
	else
		{
			verticalMVMult(in, out);
		}
}

void TSFMultiVectorOperator::verticalMVMult(const TSFVector& in, 
																						TSFVector& out) const
{
	const DenseSerialVector& xIn = TSFSerialVector::getConcrete(in);

	out.zero();
	for (int i=0; i<vectors_.length(); i++)
		{
			out.axpy(xIn[i], vectors_[i], out);
		}
}

void TSFMultiVectorOperator::horizontalMVMult(const TSFVector& in, 
																							TSFVector& out) const
{
	out = range().createMember();
	DenseSerialVector& xOut = TSFSerialVector::getConcrete(out);

	for (int i=0; i<vectors_.length(); i++)
		{
			xOut[i] = vectors_[i] * in;
		}
}

void TSFMultiVectorOperator::print(ostream& os) const 
{
	os << "<TSFMultiVectorOperator n=" << vectors_.length() << ">" << endl;

	for (int i=0; i<vectors_.length(); i++)
		{
			os << "<Vector i=" << i << ">" << endl;
			os << vectors_[i] << endl;
			os << "</Vector>" << endl;
		}
	os << "</TSFMultiVectorOperator>" << endl;
	
}
