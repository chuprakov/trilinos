#include "TSFBlockVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorSpaceBase.h"

#include "TSFProductSpace.h"
#include "TSFError.h"
#include "TSFUtils.h"

using namespace TSF;
using namespace std;

TSFBlockVector::TSFBlockVector(const TSFVectorSpace& space)
	: TSFVectorBase(space), subvectors_(space.numBlocks())
{
	for (int i=0; i<space.numBlocks(); i++)
		{
			subvectors_[i] = space.getBlock(i).createMember();
		}
}

int TSFBlockVector::numBlocks() const 
{
	return subvectors_.size();
}

void TSFBlockVector::getBlock(int i, const TSFVector& /* self */, 
																TSFVector& sub) const 
{
	sub = subvectors_[i];
}

void TSFBlockVector::setBlock(int i, const TSFVector& sub)
{
	subvectors_[i] = sub;
}

TSFReal& TSFBlockVector::setElement(int /* g */)
{
	TSFError::raise("TSFBlockVector::setElement should not be called. Set the"
									"elements in each block");
	return dummyElement_;
}

const TSFReal& TSFBlockVector::getElement(int /* g */) const 
{
	TSFError::raise("TSFBlockVector::getElement should not be called. Get the"
									"elements from each block");
	return dummyElement_; // -Wall
}

void TSFBlockVector::setElements(int /* n */, const int* /* globalIndices */,
																 const TSFReal* /* values */) 
{
	TSFError::raise("TSFBlockVector::setElements should not be called.");
}

void TSFBlockVector::getElements(int /* n */, const int* /* globalIndices */,
																 TSFReal* /* values */) const 
{
	TSFError::raise("TSFBlockVector::getElements should not be called. Get the"
									"elements from each block");
}

void TSFBlockVector::addToElements(int n, const int* /* globalIndices */,
																	 const TSFReal* /* values */) 
{	
	TSFError::raise("TSFBlockVector::addToElements should not be called. Get the"
									"elements from each block");
}

void TSFBlockVector::axpy(const TSFReal& a, const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingAxpy(a, other.getBlock(i));
		}
}

void TSFBlockVector::acceptCopyOf(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].acceptCopyOf(other);
		}
}

void TSFBlockVector::scalarMult(const TSFReal& a)
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingScalarMult(a);
		}
}

void TSFBlockVector::dotStar(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingDotStar(other.getBlock(i));
		}
}

void TSFBlockVector::dotSlash(const TSFVector& other) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].selfModifyingDotSlash(other.getBlock(i));
		}
}


TSFReal TSFBlockVector::dot(const TSFVector& other) const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].dot(other.getBlock(i));
		}
	return sum;
}

TSFReal TSFBlockVector::norm1() const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].norm1();
		}
	return sum;
}

TSFReal TSFBlockVector::normInf() const 
{
	TSFReal biggest = TSFUtils::negativeInfinity();
	for (int i=0; i<numBlocks(); i++)
		{
			TSFReal x = subvectors_[i].normInf();
			if (biggest < x) biggest = x;
		}
	return biggest;
}

TSFReal TSFBlockVector::sumElements() const
{
	TSFReal sum = 0.0;

	for (int i=0; i<numBlocks(); i++)
		{
			sum += subvectors_[i].sumElements();
		}
	return sum;
}

void TSFBlockVector::setScalar(const TSFReal& a) 
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].setScalar(a);
		}
}

TSFReal TSFBlockVector::findExtremeValue(MinOrMax type, TSFGeneralizedIndex& location, 
																				 const TSFReal& tol) const
{
	TSFReal current;
	if (type==MIN) current = TSFUtils::infinity();
	else current = TSFUtils::negativeInfinity();

	for (int i=0; i<numBlocks(); i++)
		{
			TSFGeneralizedIndex j;
			TSFReal blockVal = subvectors_[i].findExtremeValue(type, j, tol);
			if (type==MIN && blockVal < current)
				{
					location = TSFGeneralizedIndex(j, i);
					current = blockVal;
				}
			if (type==MAX && blockVal > current)
				{
					location = TSFGeneralizedIndex(j, i);
					current = blockVal;
				}
		}
	return current;
}

void TSFBlockVector::randomize(const TSFRandomNumberGenerator& r)
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].randomize(r);
		}
}

void TSFBlockVector::abs()
{
	for (int i=0; i<numBlocks(); i++)
		{
			subvectors_[i].abs();
		}
}

TSFVectorBase* TSFBlockVector::deepCopy() const 
{
	TSFVectorBase* rtn = new TSFBlockVector(space());
	for (int i=0; i<numBlocks(); i++)
		{
			rtn->setBlock(i, subvectors_[i].copy());
		}
	return rtn;
}

void TSFBlockVector::synchronizeGhostValues() const
{
	for (int i=0; i<numBlocks(); i++) subvectors_[i].synchronizeGhostValues();
}

void TSFBlockVector::invalidateGhostValues() 
{
	for (int i=0; i<numBlocks(); i++) subvectors_[i].invalidateGhostValues();
}

ostream& TSFBlockVector::print(ostream& os) const 
{
	for (int i=0; i<numBlocks(); i++) os << subvectors_[i] << endl;
	return os;
}

