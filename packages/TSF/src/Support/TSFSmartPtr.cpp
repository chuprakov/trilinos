#include "TSFSmartPtr.h"
#include "TSFError.h"

using namespace TSF;

namespace TSF
{
	void smartPtrError(const string& msg)
	{
		TSFError::raise(msg);
	}
}
