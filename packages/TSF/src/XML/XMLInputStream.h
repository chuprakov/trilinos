#ifndef XMLINPUTSTREAM_H
#define XMLINPUTSTREAM_H

#include "TSFDefs.h"

namespace TSF
{

  class XMLInputStream
    {
    public:
      /** ctor */
      XMLInputStream(){;}
      /** virtual dtor */
      inline virtual ~XMLInputStream(){;}

      /** read up to maxToRead bytes from the stream */
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead) = 0 ;

      /** identify current position */
      virtual unsigned int curPos() const ;

    };
}
#endif

