#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSource.h"
#include "TSFOperatorSourceBase.h"

using namespace TSF;
using namespace SPP;
using std::string;

RightBlockNSOperatorSource::
RightBlockNSOperatorSource(){;}

TSFLinearOperator RightBlockNSOperatorSource
::getOp() const
{ return 0; }

string RightBlockNSOperatorSource::toString() const 
{
	return "Right Block NS operator source";
}
