#include "TSFVectorSpaceBase.h"
#include "TSFError.h"
#include "TSFVectorSpace.h"

using namespace TSF;


bool TSFVectorSpaceBase::checkEquality(const TSFVectorSpaceBase* other) const
{
	/* if this method is being called, we've already tested for equality
	 * of type and dimension. The default implementation simply returns true.
	 * Some derived classes might need to make a more involved check. */

	return true;
}

void TSFVectorSpaceBase::getBlock(int i, const TSFVectorSpace& self, 
																	TSFVectorSpace& rtn) const 
{
	if (i==0) rtn = self;
	else TSFError::raise("index error in TSFVectorSpaceBase::getBlock");
}


int TSFVectorSpaceBase::typeIDCounter_ = 0 ;

int TSFVectorSpaceBase::getNewID()
{
	return typeIDCounter_++;
}

