#include "StringInputSource.h"
#include "StringInputStream.h"

using namespace TSF;


StringInputSource::StringInputSource(const string& text)
	: XMLInputSource(), text_(text)
{;}

TSFSmartPtr<XMLInputStream> StringInputSource::stream() const 
{
	return new StringInputStream(text_);
}

