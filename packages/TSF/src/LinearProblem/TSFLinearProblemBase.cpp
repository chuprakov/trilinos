#include "TSFLinearProblemBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFLinearOperatorBase.h" 
#include "TSFError.h"

using namespace TSF;

TSFVector TSFLinearProblemBase::getKnownSolution() const 
{
	TSFError::raise("getKnownSolution() called for a problem w/o "
									"known solution");
	return TSFVector(0); // -Wall
}

TSFVector TSFLinearProblemBase::getKnownSolution(const TSFVectorSpace& /* space*/) const 
{
	TSFError::raise("getKnownSolution() called for a problem w/o "
									"known solution");
	return TSFVector(0); // -Wall
}

TSFVector TSFLinearProblemBase::getRHS() const 
{
	TSFError::raise("getRHS() called for a problem w/o "
									"a specification of the vector space");
	return TSFVector(0); // -Wall
}


