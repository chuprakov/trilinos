#include "TSFUnsharedPtr.h"
#include "TSFError.h"

using namespace TSF;

namespace TSF
{
	void unsharedPtrError(const string& msg)
	{
		TSFError::raise(msg);
	}
}
