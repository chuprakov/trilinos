#include "TSFDefaultWriter.h"
#include "TSFUtils.h"
#include "TSFMPI.h"
#include "StrUtils.h"

using namespace TSF;

string& TSFDefaultWriter::header()
{
	static string rtn = "<p=";
	return rtn;
}

TSFDefaultWriter::TSFDefaultWriter()
	: os_(std::cerr)
{}

TSFDefaultWriter::TSFDefaultWriter(std::ostream& os)
	: os_(os)
{}

TSFDefaultWriter::TSFDefaultWriter(const TSFSmartPtr<ostream>& os_ptr)
	: os_ptr_(os_ptr), os_(*os_ptr_)
{}

void TSFDefaultWriter::print(const std::string& msg)
{
	string head = header() + TSFUtils::toString(TSFMPI::getRank()) + "> ";
	unsigned int maxLineSize = 78 - head.length();

	TSFArray<string> tokens = StrUtils::getTokensPlusWhitespace(msg);
	
	unsigned int lineSize = 0;
	os_ << head;
	for (int i=0; i<tokens.length(); i++)
		{
			if (lineSize+tokens[i].length() > maxLineSize)
				{
					if (StrUtils::isWhite(tokens[i])) continue;
					os_ << std::endl << head << tokens[i];
					lineSize = 0;
				}
			else
				{
					os_ << tokens[i];
				}
			lineSize += tokens[i].length();
		}
}

void TSFDefaultWriter::println(const std::string& msg)
{
	print(msg);
	os_ << std::endl;
}
