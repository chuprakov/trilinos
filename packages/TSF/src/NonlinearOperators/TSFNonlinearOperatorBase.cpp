#include "TSFError.h"
#include "TSFNonlinearOperator.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"


using namespace TSF;
using std::ostream;

TSFLinearOperator TSFNonlinearOperatorBase::derivative(const TSFVector& evalPt) const
{
	TSFError::raise("TSFNonlinearOperatorBase::derivative not available for "
									"base class");
	return TSFLinearOperator();
}
