#include "FunctionValueConvergenceTest.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpace.h"
#include "TSFVector.h"


using namespace TSF;

FunctionValueConvergenceTest::FunctionValueConvergenceTest(const TSFReal& fTol) : fTol_(fTol) 
{
	if (fTol_ < 0.0) TSFError::raise("FunctionValueConvergenceTest ctor given negative tolerance");
}

bool FunctionValueConvergenceTest::testFunctionValue(const TSFVector& f) const 
{
	return (f.norm2() < fTol_);
}
