#include "TSFVectorSpace.h"
#include "TSFVector.h"
#include "TSFEmptyVectorSpace.h"

#include "TSFError.h"
#include <typeinfo>

using namespace TSF;
using std::string;
using std::ostream;

TSFVectorSpace::TSFVectorSpace()
	: ptr_(new TSFEmptyVectorSpace())
{}

TSFVectorSpace::TSFVectorSpace(TSFVectorSpaceBase* ptr)
	: ptr_(ptr)
{;}

bool TSFVectorSpace::operator==(const TSFVectorSpace& other) const 
{
	/* first make sure our two spaces have the same type and dimension */
	if (typeid(*other.ptr_) != typeid(*ptr_) || other.dim() != dim()) return false;

	/* if they're the same type and dimension, defer to the
	 * implementations to decide if they are in fact equivalent spaces */
	return ptr_->checkEquality(other.ptr_);
}

bool TSFVectorSpace::operator!=(const TSFVectorSpace& other) const
{
	return !(operator==(other));
}

TSFVectorSpace TSFVectorSpace::deepCopy() const 
{
	return ptr_->deepCopy();
}

int TSFVectorSpace::dim() const 
{
	return ptr_->dim();
}

TSFVectorBase* TSFVectorSpace::createMember() const 
{
	return ptr_->createMember(*this);
}

bool TSFVectorSpace::contains(const TSFVector& vec) const 
{
	return (vec.space() == *this) ;
}

int TSFVectorSpace::numBlocks() const 
{
	return ptr_->numBlocks();
}

TSFVectorSpace TSFVectorSpace::getBlock(int i) const 
{
	TSFVectorSpace rtn;
	if (i<0 || i>=numBlocks()) 
		{
			TSFError::raise("index error in TSFVectorSpace::getBlock");
		}
	ptr_->getBlock(i, *this, rtn);
	return rtn;
	
}


ostream& TSFVectorSpace::print(ostream& os) const 
{
	return ptr_->print(os);
}


