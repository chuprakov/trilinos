#include "TSFUtils.h"
#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFZeroOperator.h"
#include "TSFSumOperator.h"
#include "TSFComposedOperator.h"
#include "TSFScaledOperator.h"
#include "TSFInverseOperator.h"
#include "TSFInverseAdjointOperator.h"
#include "TSFAdjointOperator.h"
#include "TSFLinearSolver.h"
#include "TSFMatrixOperator.h"
#include "TSFAdjointMatrixOperator.h"
#include "TSFOut.h"
#include "TSFUtils.h"
#include <string>


using namespace TSF;

TSFTimer TSFLinearOperator::opTimer_("Linear operators");

TSFLinearOperator::TSFLinearOperator()
  :ptr_(0), hasTranspose_(false)
{;}

TSFLinearOperator::TSFLinearOperator(TSFLinearOperatorBase* ptr)
	: ptr_(ptr), hasTranspose_(false)
{;}

const TSFVectorSpace& TSFLinearOperator::range() const 
{
	return ptr_->range();
}

const TSFVectorSpace& TSFLinearOperator::domain() const 
{
	return ptr_->domain();
}

bool TSFLinearOperator::isZeroOperator() const 
{
	return ptr_->isZeroOperator();
}

int TSFLinearOperator::numBlockRows() const 
{
	return ptr_->numBlockRows();
}

int TSFLinearOperator::numBlockCols() const 
{
	return ptr_->numBlockCols();
}

bool TSFLinearOperator::isBlockOperator() const
{
  return ptr_->isBlockOperator();
}

TSFLinearOperator TSFLinearOperator::getBlock(int i, int j) const
{
	if (ptr_->isBlockOperator()) 
		{
			TSFLinearOperator rtn;
			ptr_->getBlock(i, j, rtn);
			return rtn;
		}
	else 
		{
			if (i==0 && j==0)
				{
					return *this;
				}
			else
				{
					TSFError::raise("TSFLinearOperator::getBlock index out of range");
					return *this;
				}
		}
}

TSFLinearOperator& TSFLinearOperator::setBlock(int i, int j,
																							 const TSFLinearOperator& sub)
{
	ptr_->setBlock(i, j, sub);
	return *this;
}



void TSFLinearOperator::getRow(const int row, TSFArray<int>& indices,
                               TSFArray<TSFReal>& values) const
{
  if (ptr_->isBlockOperator())
    {
      /* find col block in which row exists and row within the block */
      int K = 0;
      int blockRow = 0;
      int rowInBlock = 0;
      for (int i = 0; i < numBlockRows(); i++)
        {
          int numR = getBlock(i, 0).range().dim();
          K += numR;
          if (row < K)
            {
              blockRow = i;
              rowInBlock = row - (K - numR);
              break;
            }
        }
/*       cerr << "In TSFLinearOp: row = " << row << " rowInBlock = "  */
/*            << rowInBlock << endl; */

      /* get the row elements for each block in the row.  */
      int offset = 0;
      TSFArray<int> localInd;
      TSFArray<TSFReal> localVal;
      for (int i = 0; i < numBlockCols(); i++)
        {
          localInd.resize(0);
          localVal.resize(0);
          TSFLinearOperator bl;
          ptr_->getBlock(blockRow, i, bl);
          bl.getRow(rowInBlock,localInd, localVal); 
          for (int j = 0; j < localInd.size(); j++)
            {
              indices.append(localInd[j] + offset);
              values.append(localVal[j]);
            }
          offset += getBlock(blockRow, i).domain().dim();
        }
    }
  else
    {
      ptr_->getRow(row, indices, values);
    }
} 

TSFLinearOperator& TSFLinearOperator::getTranspose() 
{
/*   cerr << "Getting transpose for \n"; */
/*   describe(); */
  if (!hasTranspose_)
    {
      //cerr << "building trp\n";
          opTrp_ = ptr_->getTranspose();
          opTrp_->setTranspose(this);
          hasTranspose_ = true;
    }
  return *opTrp_;
}

void TSFLinearOperator::setTranspose(TSFLinearOperator* op)
{
  opTrp_ = op;
}


void TSFLinearOperator::apply(const TSFVector& arg,
															TSFVector& out) const 
{
	TSFTimeMonitor t(opTimer_);
	if (out.isNull()) out = range().createMember();
	else (out.zero());
	ptr_->apply(arg, out);
}

void TSFLinearOperator::applyInverse(const TSFVector& arg,
																		 TSFVector& out) const 
{
	TSFTimeMonitor t(opTimer_);
	if (out.isNull()) out = domain().createMember();
	ptr_->applyInverse(arg, out);
}

void TSFLinearOperator::applyInverse(const TSFLinearSolver& solver,
																		 const TSFVector& arg,
																		 TSFVector& out) const 
{
	TSFTimeMonitor t(opTimer_);
	if (out.isNull()) out = domain().createMember();
	solver.solve(*this, arg, out);
}

void TSFLinearOperator::applyAdjoint(const TSFVector& arg,
																		 TSFVector& out) const 
{
	TSFTimeMonitor t(opTimer_);
	if (out.isNull()) out = domain().createMember();
	ptr_->applyAdjoint(arg, out);
}

void TSFLinearOperator::applyInverseAdjoint(const TSFVector& arg,
																						TSFVector& out) const 
{
	TSFTimeMonitor t(opTimer_);
	if (out.isNull()) out = range().createMember();
	ptr_->applyInverseAdjoint(arg, out);
}

void TSFLinearOperator::applyInverseAdjoint(const TSFLinearSolver& solver,
																						const TSFVector& arg,
																						TSFVector& out) const 
{
	if (out.isNull()) out = range().createMember();
	solver.solve(adjoint(), arg, out);
}

TSFVector TSFLinearOperator::operator*(const TSFVector& x) const 
{
	TSFVector rtn = range().createMember();
/*     cerr << "In TSFLinearOperator:*\n"; */
/*     rtn.describe(); */
	apply(x, rtn);
	return rtn;
}

TSFLinearOperator TSFLinearOperator::inverse() const 
{
	TSFLinearOperator rtn;
	ptr_->getInverse(*this, rtn);
	return rtn;
}

TSFLinearOperator TSFLinearOperator::inverse(const TSFLinearSolver& solver) const 
{
	TSFLinearOperator rtn;
	ptr_->getInverse(solver, *this, rtn);
	return rtn;
}

TSFLinearOperator TSFLinearOperator::adjoint() const 
{
	TSFLinearOperator rtn;
	ptr_->getAdjoint(*this, rtn);
	return rtn;
}

TSFLinearOperator TSFLinearOperator::inverseAdjoint() const 
{
	TSFLinearOperator rtn;
	ptr_->getInverseAdjoint(*this, rtn);
	return rtn;
}

TSFLinearOperator TSFLinearOperator::inverseAdjoint(const TSFLinearSolver& solver) const 
{
	TSFLinearOperator rtn;
	ptr_->getInverseAdjoint(solver, *this, rtn);
	return rtn;
}


TSFLinearOperator TSFLinearOperator::operator+(const TSFLinearOperator& op) const 
{
	/* do consistency checking on spaces */
	if (range() != op.range())
		{
			TSFError::raise("range mismatch in TSFLinearOperator::operator+");
		}
	if (domain() != op.domain())
		{
			TSFError::raise("domain mismatch in TSFLinearOperator::operator+");
		}

	/* handle the special case of addition with zero */
	if (isZeroOperator())
		{
			return op;
		}
	else if (op.isZeroOperator())
		{
			return *this;
		}
	/* in general case, create a sum operator */
	return new TSFSumOperator(*this, op);
}

TSFLinearOperator TSFLinearOperator::operator-(const TSFLinearOperator& op) const 
{
	/* do consistency checking on spaces */
	if (range() != op.range())
		{
			TSFError::raise("range mismatch in TSFLinearOperator::operator-");
		}
	if (domain() != op.domain())
		{
			TSFError::raise("domain mismatch in TSFLinearOperator::operator");
		}

	/* handle the special case of subtraction with zero */
	if (isZeroOperator())
		{
			return -op;
		}
	else if (op.isZeroOperator())
		{
			return *this;
		}
	/* in general case, create a sum operator with a boolean to 
	 * indicate subtraction*/
	return new TSFSumOperator(*this, op, true);
}

TSFLinearOperator TSFLinearOperator::operator*(const TSFLinearOperator& op) const 
{
	/* check for consistency between the range of the righthand operator
	* and the domain of the lefthand operator */
	if (op.range() != domain())
		{
			TSFError::raise("domain-range mismatch in TSFLinearOperator::operator*");
		}

	/* check for composition with a zero operator. If one of the operands is
	* zero, create a new zero operator with the appropriate domain and range */
	if (op.isZeroOperator())
		{
			return new TSFZeroOperator(op.domain(), range());
		}
	else if (isZeroOperator())
		{
			return new TSFZeroOperator(op.domain(), range());
		}

	/* general case, return a composed operator */
	return new TSFComposedOperator(*this, op);
}

TSFLinearOperator TSFLinearOperator::operator*(const TSFReal& scale) const 
{
	/* check for special cases op=0, scale=0, and scale=1 */
	if (isZeroOperator())
		{
			return *this;
		}
	else if (TSFUtils::chop(scale)==0)
		{
			return new TSFZeroOperator(domain(), range());
		}
	else if (TSFUtils::chop(scale-1.0)==0)
		{
			return *this;
		}	
	/* general case, return a scaled operator */
	return new TSFScaledOperator(*this, scale);
}

TSFLinearOperator TSFLinearOperator::operator-() const 
{
	return new TSFScaledOperator(*this, -1.0);
}

bool TSFLinearOperator::isMatrixOperator() const 
{
	return ptr_->isMatrixOperator();
}

const TSFSmartPtr<const TSFMatrixOperator> TSFLinearOperator::getMatrix() const 
{
	const TSFSmartPtr<const TSFMatrixOperator>
		mat_op = ptr_->getMatrix();
	if (mat_op.isNull())
		{
			TSFError::raise("TSFLinearOperator::getMatrix() called for matrix-free operator");
		}
	return mat_op;
}

/* TSFSmartPtr<TSFMatrixOperator> TSFLinearOperator::getMatrixC()   */
/* { */
/* 	TSFSmartPtr<TSFMatrixOperator> */
/* 		mat_op = ptr_->getMatrixC(); */
/* 	if (mat_op.isNull()) */
/* 		{ */
/* 			TSFError::raise("TSFLinearOperator::getMatrixC() called for matrix-free operator"); */
/* 		} */
/* 	return mat_op; */
/* } */

string TSFLinearOperator::toString() const 
{
  TSFOStringStream ost;
            ost.clear();
	ptr_->print(ost);
    ost << ends;
    
	return ost.str();
}


/* Describe the operator in terms of Blocks  */
void TSFLinearOperator::describe() const
{
  describe(0);
}

/* Indent to get nice spacing */
void TSFLinearOperator::describe(const int&depth) const
{
  string spaces = "";
  for (int i = 0; i < depth; i++)
    {
      spaces.append("  ");
    }

  if (ptr_->isBlockOperator())
    {
      TSFOut::println(spaces + "Block Operator with " 
                      + TSFUtils::toString(numBlockRows())
                      + " rows and  "  + TSFUtils::toString(numBlockCols())   
                      + " columns");
/*       cerr << spaces << "Block Operator with " << numBlockRows() */
/*            << " rows and  " << numBlockCols()   */
/*            << " columns" << endl; */
      for (int i = 0; i < numBlockRows(); i++)
        {
          for (int j = 0; j < numBlockCols(); j++)
            {
              TSFOut::println(spaces + "  Block " + TSFUtils::toString(i) + ", " 
                              + TSFUtils::toString(j) + " is");
              getBlock(i, j).describe(depth+1);
              TSFOut::println(spaces + "  Block " + TSFUtils::toString(i) + ", " 
                              + TSFUtils::toString(j) + " is");

/*               cerr << spaces << "  Block " << i << ", "  */
/*               << j << " is\n "; */
              getBlock(i, j).describe(depth+1);
            }
        }
    }
  else
    {
      string typeNam = ptr_->typeName();
      //string typeNam = "junk";
      TSFOut::println(spaces + "  " + typeNam + " of dimension (rows) " 
           + TSFUtils::toString(range().dim())  + " by  (cols) " 
           + TSFUtils::toString(domain().dim()));

/*       cerr << spaces << typeNam <<" of dimension "  */
/*            << domain().dim() << " by  "  */
/*            << range().dim() << endl;  */
    }
}
  




namespace TSF
{
	TSFLinearOperator operator*(const TSFReal& scale, 
															const TSFLinearOperator& op)
	{
		return op*scale;
	}

	ostream& operator<<(ostream& os, const TSFLinearOperator& op)
	{
		op.ptr_->print(os);
		return os;
	}
}


 
