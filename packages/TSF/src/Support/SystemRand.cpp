#include "SystemRand.h"

#include <stdlib.h>

using namespace TSF;

SystemRand::SystemRand()
{
	srand(TSFRandomNumberGeneratorBase::microsecondSeed());
}

SystemRand::SystemRand(unsigned int seed)
{
	srand(seed);
}

void SystemRand::generateRandomNumbers(int n, TSFReal* x) const
{
	for (int i=0; i<n; i++)
		{
			*x++ = ((double) rand()) / ((double) RAND_MAX);
		}
}



