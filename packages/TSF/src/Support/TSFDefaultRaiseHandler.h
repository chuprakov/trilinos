#ifndef TSFDEFAULTRAISEHANDLER_H
#define TSFDEFAULTRAISEHANDLER_H

#include "TSFConfig.h"
#include "TSFRaiseHandlerBase.h"
#include <stdexcept>

namespace TSF
{
	using std::string;

	/** \ingroup ErrorHandling
	 * The default raise handler throws a std::runtime_error() exception.
	 */
	class TSFDefaultRaiseHandler : public TSFRaiseHandlerBase
		{
		public:
			/** empty ctor */
			TSFDefaultRaiseHandler(){;}

			/** TUVD */
			virtual ~TSFDefaultRaiseHandler(){;}

			/** When TSFError::raise() is called, throw a std::runtime_error 
			 * exception. */
			virtual void handleRaise(const char* msg);

		private:
		};
}

#endif
