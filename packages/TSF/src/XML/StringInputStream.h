#ifndef STRINGINPUTSTREAM_H
#define STRINGINPUTSTREAM_H

#include "TSFConfig.h"

#include "XMLInputStream.h"
#include <cstdio>
#include <string>

namespace TSF
{
	using std::string;

	/** \ingroup XML 
	 * Input stream for reading an entire document from a string. This
	 * is a low-level object and should not be needed at the user level.
	 * FileInputSource is the user-level object.
	 */
	
	class StringInputStream : public XMLInputStream
		{
		public:
			/** construct with the string from which data will be read */
			StringInputStream(const string& text) 
				: XMLInputStream(), text_(text), pos_(0) {;}
			
			/** virtual dtor */
			virtual ~StringInputStream() {;}
			
			/** read up to maxToRead bytes */
			virtual unsigned int readBytes(unsigned char* const toFill, 
																		 const unsigned int maxToRead);
			
		private:
			string text_;
			unsigned int pos_;
		};

}
#endif

