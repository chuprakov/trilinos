#ifndef TSFDEFAULTWRITER_H
#define TSFDEFAULTWRITER_H

#include "TSFConfig.h"
#include "TSFWriterBase.h"
#include <iostream>

namespace TSF
{
	using std::string;
	using std::ostream;

	/** \ingroup IO
	 * 
	 */
	class TSFDefaultWriter : public TSFWriterBase
		{
		public:
			/** */
			TSFDefaultWriter();

			/** */
			TSFDefaultWriter(ostream& os);

			/** */
			virtual ~TSFDefaultWriter(){;}

			/** */
			virtual void print(const string& msg);

			/** */
			virtual void println(const string& msg);

		private:
			std::ostream& os_;
			static string& header();
		};
}
#endif
