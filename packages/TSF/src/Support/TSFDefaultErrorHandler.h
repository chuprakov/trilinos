#ifndef TSFDEFAULTERRORHANDLER_H
#define TSFDEFAULTERRORHANDLER_H

#include "TSFConfig.h"
#include <stdexcept>

namespace TSF
{
	using std::string;
	
	/** \ingroup ErrorHandling
	 *
	 */
	class TSFDefaultErrorHandler
		{
		public:
			/** empty ctor */
			TSFDefaultErrorHandler(){;}

			/** TUVD */
			virtual ~TSFDefaultErrorHandler(){;}

			/** */
			virtual void handleError(const char* msg);

		private:
		
		};
}
#endif
