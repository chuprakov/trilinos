#include "TSFDeferredLinearCombination.h"

#include "TSFUtils.h"

using namespace TSF;

TSFDeferredLinearCombination::TSFDeferredLinearCombination(const TSFReal& coeff, 
																					 const TSFVector& v) 
	: coeffs_(1), vectors_(1)
{
	coeffs_[0] = coeff;
	vectors_[0] = v;
}

void TSFDeferredLinearCombination::add(const TSFReal& coeff, 
																			 const TSFVector& v)
{
	coeffs_.append(coeff);
	vectors_.append(v);
}

TSFDeferredLinearCombination TSFDeferredLinearCombination::operator-() const 
{
	TSFDeferredLinearCombination rtn = *this;
	for ( int i=0; i<rtn.coeffs_.size(); i++)
		{
			rtn.coeffs_[i] *= -1.0;
		}
	return rtn;
}

TSFDeferredLinearCombination TSFDeferredLinearCombination::operator+(const TSFDeferredLinearCombination& op) const 
{
	TSFDeferredLinearCombination rtn = *this;
	
	for ( int i=0; i<op.coeffs_.size(); i++)
		{
			rtn.add(op.coeffs_[i], op.vectors_[i]);
		}
	
	return rtn;
}

TSFDeferredLinearCombination TSFDeferredLinearCombination::operator-(const TSFDeferredLinearCombination& op) const 
{
	TSFDeferredLinearCombination rtn = *this;
	
	for ( int i=0; i<op.coeffs_.size(); i++)
		{
			rtn.add(-op.coeffs_[i], op.vectors_[i]);
		}
	
	return rtn;
}


TSFDeferredLinearCombination TSFDeferredLinearCombination::operator+(const TSFVector& v) const 
{
	TSFDeferredLinearCombination rtn = *this;
	
	rtn.add(1.0, v);

	return rtn;
}

TSFDeferredLinearCombination TSFDeferredLinearCombination::operator-(const TSFVector& v) const 
{
	TSFDeferredLinearCombination rtn = *this;
	
	rtn.add(-1.0, v);

	return rtn;
}

TSFDeferredLinearCombination TSFDeferredLinearCombination::operator*(const TSFReal& a) const 
{
	TSFDeferredLinearCombination rtn = *this;

	for ( int i=0; i<rtn.coeffs_.size(); i++)
		{
			rtn.coeffs_[i] *= a;
		}
	return rtn;
}

TSFVector TSFDeferredLinearCombination::evaluateIntoNewVector() const 
{
	TSFVector rtn = vectors_[0].copy();
	evaluateIntoExistingVector(rtn);
	return rtn;
}

void TSFDeferredLinearCombination::evaluateIntoExistingVector(TSFVector& result) const
{
	/* We have to be a bit careful here. The target vector may also be one of the
	 * terms in the sum. In that case, we must be sure to start the accumulation
	 * with that term. */

	TSFArray<bool> LHSVectorFlag = markLHSVectors(result);

	if (LHSVectorFlag.size()==0)
		{
			TSFReal aStart = coeffs_[0];
			if (TSFUtils::chop(aStart-1.0)!=0) 
				{
					result.scalarMult(aStart, vectors_[0]);
				}
			else
				{
					result = vectors_[0].copy();
				}
		}
	else
		{
			int pos = 0;
			TSFReal coeff = 0.0;
			for ( int i=0; i<coeffs_.size(); i++)
				{
					if (!LHSVectorFlag[i]) continue;
					pos=i;
					coeff += coeffs_[i];
				}
			if (TSFUtils::chop(coeff-1.0)!=0) 
				{
					result.scalarMult(coeff, vectors_[pos]);
				}
			else
				{
					result = vectors_[pos].copy();
				}
		}

	/* accumulate the remaining terms into the sum */
	for ( int i=0; i<coeffs_.size(); i++)
		{
			if (LHSVectorFlag[i]) continue;
			result.selfModifyingAxpy(coeffs_[i], vectors_[i]);
		}
}

TSFReal TSFDeferredLinearCombination::operator*(const TSFDeferredLinearCombination& b) const 
{
	if (coeffs_.size()==1)
		{
			return coeffs_[0] * (vectors_[0] * b);
		}
	else
		{
			return evaluateIntoNewVector() * b.evaluateIntoNewVector();
		}
}

TSFReal TSFDeferredLinearCombination::operator*(const TSFVector& b) const 
{
	if (coeffs_.size()==1)
		{
			return coeffs_[0] * (vectors_[0] * b);
		}
	else
		{
			return evaluateIntoNewVector() * b;
		}
}

TSFArray<bool> TSFDeferredLinearCombination::markLHSVectors(const TSFVector& lhs) const 
{
	TSFArray<bool> rtn(vectors_.size());
	for ( int i=0; i<vectors_.size(); i++)
		{
			rtn[i] = lhs.isIdenticalTo(vectors_[i]);
		}
	return rtn;
}

TSFReal TSFDeferredLinearCombination::norm2() const 
{
	if (coeffs_.size()==1)
		{
			return coeffs_[0] * vectors_[0].norm2();
		}
	else
		{
			return evaluateIntoNewVector().norm2();
		}
}

TSFReal TSFDeferredLinearCombination::norm1() const 
{
	if (coeffs_.size()==1)
		{
			return coeffs_[0] * vectors_[0].norm1();
		}
	else
		{
			return evaluateIntoNewVector().norm1();
		}
}

TSFReal TSFDeferredLinearCombination::normInf() const 
{
	if (coeffs_.size()==1)
		{
			return coeffs_[0] * vectors_[0].normInf();
		}
	else
		{
			return evaluateIntoNewVector().normInf();
		}
}


namespace TSF
{
	TSFDeferredLinearCombination operator+(const TSFVector& a, 
																				 const TSFVector& b)
	{
		TSFDeferredLinearCombination rtn(1.0, a);
		return rtn + b;
	}

	TSFDeferredLinearCombination operator-(const TSFVector& a, 
																				 const TSFVector& b)
	{
		TSFDeferredLinearCombination rtn(1.0, a);
		return rtn - b;
	}

	
}

