#include "TSFSerialVector.h"
#include "TSFError.h"

using namespace TSF;

TSFSerialVector::TSFSerialVector(const TSFVectorSpace& space)
	: TSFInCoreVector(space), x_(space.dim())
{}

void TSFSerialVector::setElements(int n, const int* globalIndices,
																	const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			x_[globalIndices[i]] = values[i];
		}
}

void TSFSerialVector::getElements(int n, const int* globalIndices,
																	TSFReal* values) const 
{
	for (int i=0; i<n; i++)
		{
			values[i] = x_[globalIndices[i]];
		}
}

void TSFSerialVector::addToElements(int n, const int* globalIndices,
																		const TSFReal* values) 
{
	for (int i=0; i<n; i++)
		{
			x_[globalIndices[i]] += values[i];
		}
}

TSFVectorBase* TSFSerialVector::deepCopy() const 
{
	TSFVectorBase* rtn = new TSFSerialVector(*this);
	if (rtn==0) TSFError::raise("TSFSerialVector::deepCopy()");
	return rtn;
}

ostream& TSFSerialVector::print(ostream& os) const 
{
	os << x_ ;
	return os;
}

const DenseSerialVector& TSFSerialVector::getConcrete(const TSFVector& x)
{
	const TSFSerialVector* v = dynamic_cast<const TSFSerialVector*>(x.ptr());
	if (v==0) TSFError::raise("bad cast in TSFSerialVector::getConcrete");
	return v->x_;
}

DenseSerialVector& TSFSerialVector::getConcrete(TSFVector& x)
{
	TSFSerialVector* v = dynamic_cast<TSFSerialVector*>(x.ptr());
	if (v==0) TSFError::raise("bad cast in TSFSerialVector::getConcrete");
	return v->x_;
}










