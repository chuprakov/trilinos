#include "TSFDefaultWriter.h"
#include "TSFUtils.h"
#include "TSFMPI.h"

using namespace TSF;

string& TSFDefaultWriter::header()
{
	static string rtn = "<TSF p=";
	return rtn;
}

TSFDefaultWriter::TSFDefaultWriter()
	: os_(std::cerr)
{}

TSFDefaultWriter::TSFDefaultWriter(std::ostream& os)
	: os_(os)
{}

void TSFDefaultWriter::print(const std::string& msg)
{
	string out = header() + TSFUtils::toString(TSFMPI::getRank()) + "> " + msg;
	os_ << out;
}

void TSFDefaultWriter::println(const std::string& msg)
{
	print(msg);
	os_ << std::endl;
}
