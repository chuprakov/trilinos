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

void DenseSerialVectorSpace::print(ostream& os) const 
{
	os << "DenseSerialVectorSpace[" << dim() << "]";
}


