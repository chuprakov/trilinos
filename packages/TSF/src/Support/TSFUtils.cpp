#include "TSFUtils.h"
#include "TSFVersion.h"
#include "TSFOut.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

using namespace TSF;


TSFReal TSFUtils::chopVal_ = 1.0e-16;

void TSFUtils::aboutBuild()
{
	TSF::showVersion();
}

TSFReal TSFUtils::chop(const TSFReal& x) 
{
	if (fabs(x) < chopVal_) return 0;
	return x;
}

string TSFUtils::toString(const int& x)
{
	char s[100];
	sprintf(s, "%d", x);
	return string(s);
}

string TSFUtils::toString(const TSFReal& x)
{
	char s[100];
	sprintf(s, "%g", x);
	return string(s);
}



