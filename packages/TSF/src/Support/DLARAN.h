#ifndef DLARAN_H
#define DLARAN_H

#include "TSFConfig.h"
#include "TSFRandomNumberGenerator.h"

namespace TSF
{
	/**\ingroup Support 
	 * DLARAN generates a vector of uniform deviates
	 * using the 48-bit algorithm from Lapack's
	 * dlaran() generator. By default, the current millisecond count is used
	 * as a seed. Alternatively, a DLARAN object can be constructed with 
	 * an specified integer seed.
	 * 
	 * See the Lapack documentation for further information.
	 * on the algorithm.
	 * @author Kevin Long
	 */
	class DLARAN : public TSFRandomNumberGeneratorBase
		{
		public:
			/** Empty ctor initializes seed to current millisecond count */
			DLARAN();
			/** Initialize to an integer seed */
			DLARAN(unsigned int seed);
  
			/** the usual virtual dtor */
			virtual ~DLARAN() {;}
  
			/** generate a random uniform unit deviate */
			virtual void generateRandomNumbers(int n, TSFReal* x) const ;

		private:
			/** create the 48-bit seed from the initial seed */
			void init(unsigned int seed);

			/** get the next number from the sequence */
			double getNext() const ;

			/* 48-bit seed */
			mutable int seed48_[4];
		};


}

#endif 
