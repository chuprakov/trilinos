#include "StepConvergenceTest.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"

#include "TSFError.h"

using namespace TSF;

StepConvergenceTest::StepConvergenceTest(const TSFReal& stepTol) 
	: stepTol_(stepTol) 
{
	if (stepTol_ < 0.0) TSFError::raise("StepConvergenceTest ctor given negative tolerance");
}

bool StepConvergenceTest::testStep(const TSFVector& deltaX) const 
{
	return (deltaX.norm2() < stepTol_);
}
