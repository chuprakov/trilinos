#ifndef TSFPARAMETERREADERBASE_H
#define TSFPARAMETERREADERBASE_H

#include "TSFConfig.h"
#include <string>
#include "TSFArray.h"

namespace TSF 
{

	using namespace std;

	class TSFParameterReaderBase
		{
		public:
			/** empty ctor */
			TSFParameterReaderBase();

			/** virtual dtor */
			virtual ~TSFParameterReaderBase();

			/** The read() method gets a parameter list. How this is done, for instance,
			 * how a text file is parsed into a parameter list, will depend on the
			 * implementation of the derived reader type. */
			virtual TSFParameterList read() const = 0 ;

		private:
		};

}

#endif
