#include "DenseSerialVectorSpace.h"

#include "TSFUtils.h"
#include "TSFSerialVector.h"


using namespace TSF;

TSFVectorSpaceBase* DenseSerialVectorSpace::deepCopy() const 
{
	return new DenseSerialVectorSpace(*this);
}

TSFVectorBase* DenseSerialVectorSpace::createMember(const TSFVectorSpace& space) const 
{
	return new TSFSerialVector(space);
}

ostream& DenseSerialVectorSpace::print(ostream& os) const 
{
	return os << "DenseSerialVectorSpace[" << dim() << "]";
}


