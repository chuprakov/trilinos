#include "TSFRandomNumberGenerator.h"

#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

using namespace TSF;

TSFRandomNumberGenerator::TSFRandomNumberGenerator(TSFRandomNumberGeneratorBase* ptr)
	: ptr_(ptr)
{;}

unsigned int TSFRandomNumberGeneratorBase::microsecondSeed()
{

	struct timeval tp;
  struct timezone tzp;
  
  gettimeofday(&tp, &tzp);

  unsigned int rtn = tp.tv_usec;
  
  return rtn;

}

