#include "TSFDefaultRaiseHandler.h"

using namespace TSF;

void TSFDefaultRaiseHandler::handleRaise(const char* msg)
{
	throw std::runtime_error(msg);
}


