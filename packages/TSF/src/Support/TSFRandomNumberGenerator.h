#ifndef TSFRANDOMNUMBERGENERATOR_H
#define TSFRANDOMNUMBERGENERATOR_H

#include "TSFConfig.h"
#include "TSFSmartPtr.h"
#include <string>
#include <stdexcept>

namespace TSF
{
	using std::string;

	/** \ingroup Support
	 * 
	 */
	class TSFRandomNumberGeneratorBase
		{
		public:
			/** Empty ctor */
			TSFRandomNumberGeneratorBase(){;}

			/** Virtual dtor */
			virtual ~TSFRandomNumberGeneratorBase(){;}

			/** generate a vector of random numbers */
			virtual void generateRandomNumbers(int n, TSFReal* x) const = 0 ;

		protected:
			static unsigned int microsecondSeed() ;
		private:
		};

	/** \ingroup Support
	 * 
	 */
	class TSFRandomNumberGenerator
		{
		public:
			/** Construct with a pointer to a concrete type */
			TSFRandomNumberGenerator(TSFRandomNumberGeneratorBase* ptr);

			/** generate a vector of random numbers */
			void generateRandomNumbers(int n, TSFReal* x) const 
				{ptr_->generateRandomNumbers(n, x);}

		private:
			TSFSmartPtr<TSFRandomNumberGeneratorBase> ptr_;
		};
}
#endif
