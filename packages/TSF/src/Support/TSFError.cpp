#include "TSFError.h"
#include "TSFDefaultRaiseHandler.h"

#include "TSFArray.h"
#include <stdexcept>
#include "TSFOut.h"



using namespace TSF;


using std::runtime_error;

/* initialize the static raise handler object to the default handler. This
 * can be changed later with a call to setRaiseHandler() */

TSFSmartPtr<TSFRaiseHandlerBase> TSFError::handler_ = new TSFDefaultRaiseHandler();


void TSFError::raise(const std::string& msg)
{
	TSFArray<int> blah;
	TSFArray<bool> blap;
	handler_->handleRaise(msg.c_str());
}

void TSFError::trace(const std::exception& e, const std::string& where)
{
	std::string msg = std::string(e.what()) + " at " + where;

	TSFOut::println(msg);

	throw runtime_error(msg);
}
