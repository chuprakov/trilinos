#include "TSFDefaultErrorHandler.h"

void TSFDefaultErrorHandler::handleError(const char* msg)
{
	std::runtime_error.raise(msg);
}
