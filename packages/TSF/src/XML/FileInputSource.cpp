#include "FileInputSource.h"
#include "FileInputStream.h"

using namespace TSF;


FileInputSource::FileInputSource(const string& filename)
	: XMLInputSource(), filename_(filename)
{;}

TSFSmartPtr<XMLInputStream> FileInputSource::stream() const 
{
	return new FileInputStream(filename_);
}

