#include "TSFMultiVectorOperator.h"

#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFUtils.h"
#include "DenseSerialVectorSpace.h"
#include "TSFSerialVector.h"
#include "TSFDeferredLinearCombination.h"
#include "TSFLinearOperator.h"

using namespace TSF;

TSFMultiVectorOperator::TSFMultiVectorOperator()
	:TSFLinearOperatorBase(TSFVectorSpace(),TSFVectorSpace()),
	 isVertical_(false),
  vectors_(0),
  haveTranspose_(false)
{}

TSFMultiVectorOperator::TSFMultiVectorOperator(const TSFVectorSpace& space,
																							 int numVectors,
																							 bool isVertical)
	:TSFLinearOperatorBase(TSFVectorSpace(),TSFVectorSpace()),
  isVertical_(isVertical), vectors_(numVectors), haveTranspose_(false)

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

/* TSFMultiVectorOperator::TSFMultiVectorOperator(const TSFVectorSpace& space,  */
/*                                                bool isVertical, */
/*                                                const TSFArray<TSFVector>& vectors) */
/* 	:TSFLinearOperatorBase(TSFVectorSpace(),TSFVectorSpace()), */
/* 	 isVertical_(isVertical), vectors_(vectors) */

/* { */
/*   int numVectors = vectors_.length(); */
/* 	if (isVertical_) */
/* 		{ */
/* 			range_ = space; */
/* 			domain_ = new DenseSerialVectorSpace(numVectors); */
/* 		} */
/* 	else */
/* 		{ */
/* 			domain_ = space; */
/* 			range_ = new DenseSerialVectorSpace(numVectors); */
/* 		} */
/* } */


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


TSFLinearOperator* TSFMultiVectorOperator::getTranspose() 
{
  if (haveTranspose_) 
    {
      return &opTrp_;
    }
  haveTranspose_ = true;
  if (isVertical_)
    {
      opTrp_ = new TSFMultiVectorOperator(range_, vectors_.size(),
                                     !isVertical_);
    }
  else
    {
      opTrp_ = new TSFMultiVectorOperator(domain_, vectors_.size(),
                                     !isVertical_); 
    }

  const TSFSmartPtr<TSFLinearOperatorBase>& base = opTrp_.getPtr();
  const TSFLinearOperatorBase* eBase = &(*base);
  const TSFMultiVectorOperator* mvop 
    = dynamic_cast<const TSFMultiVectorOperator*>(eBase);
      
  for (int i = 0; i < vectors_.size(); i++)
    {
      mvop->getVector(0) = vectors_[i].copy();
    }
  return &opTrp_;
}


void TSFMultiVectorOperator::verticalMVMult(const TSFVector& in, 
                                            TSFVector& out) const
{
  try
    {
      const DenseSerialVector& xIn = TSFSerialVector::getConcrete(in);
      
      out.zero();
      for (int i=0; i<vectors_.length(); i++)
		{
          out.axpy(xIn[i], vectors_[i], out);
		}
    }
  catch(exception& e)
    {
      TSFError::trace(e, "in TSFMultiVectorOperator::verticalMVMult()");
    }
}

void TSFMultiVectorOperator::horizontalMVMult(const TSFVector& in, 
                                              TSFVector& out) const
{
  try
    {
      //out = range().createMember();
      DenseSerialVector& xOut = TSFSerialVector::getConcrete(out);

	for (int i=0; i<vectors_.length(); i++)
		{
			xOut[i] = vectors_[i] * in;
		}
    }
  catch(exception& e)
    {
      TSFError::trace(e, "in TSFMultiVectorOperator::horizontalMVMult()");
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


void TSFMultiVectorOperator::getRow(int row, TSFArray<int>& indices, 
                  TSFArray<TSFReal>& values) const
{
  if (isVertical_)
    {
      for (int i = 0; i < vectors_.size(); i++)
        {
          indices.append(i);
          values.append(vectors_[i][row]);
        }
    }
  else
    {
      for (int i = 0; i < vectors_[0].space().dim(); i++)
        {
          indices.append(i);
          values.append(vectors_[row][i]);
        }
    }
}
