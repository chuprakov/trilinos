#ifndef SYSTEMRAND_H
#define SYSTEMRAND_H

#include "TSFConfig.h"
#include "TSFRandomNumberGenerator.h"

namespace TSF
{
	/**\ingroup Support 
	 * Use the system's rand() function to fill a random vector. 
	 *
	 * @author Kevin Long
	 */
	class SystemRand : public TSFRandomNumberGeneratorBase
		{
		public:
			/** Empty ctor uses system's default seed */
			SystemRand();
			/** Initialize to an integer seed */
			SystemRand(unsigned int seed);
  
			/** the usual virtual dtor */
			virtual ~SystemRand() {;}
  
			/** generate a random uniform unit deviate */
			virtual void generateRandomNumbers(int n, TSFReal* x) const ;

		private:
		};


}

#endif 
