#ifndef TSFWRITERBASE_H
#define TSFWRITERBASE_H

#include "TSFConfig.h"
#include <string>
#include <stdexcept>

namespace TSF
{
	using std::string;

	/** \ingroup IO
	 * 
	 */
	class TSFWriterBase
		{
		public:
			/** Empty ctor */
			TSFWriterBase(){;}

			/** TUVD */
			virtual ~TSFWriterBase(){;}

			/** print a string */
			virtual void print(const string& msg) = 0 ;

			/** print a string followed by a newline */
			virtual void println(const string& msg) = 0 ;

		private:
		};
}
#endif
