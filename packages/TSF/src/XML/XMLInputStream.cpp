#include "XMLInputStream.h"
#include "TSFError.h"

using namespace TSF;


unsigned int XMLInputStream::curPos() const 
{
	TSFError::raise("XMLInputStream::curPos() should never be called. It is"
							 "there only for compatibility with Xerces");
	return 0;
}
