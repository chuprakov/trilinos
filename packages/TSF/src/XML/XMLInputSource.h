#ifndef XMLINPUTSOURCE_H
#define XMLINPUTSOURCE_H

#include "TSFConfig.h"
#include "XMLObject.h"
#include "XMLInputStream.h"

namespace TSF
{
	/** \ingroup XML
	 * XMLInputSource represents a source of XML input that can be parsed
	 * to produce an XMLObject. The source might be a file, a socket, a
	 * string. The XMLObject is created with a call to the getObject() method.
	 *
	 * The source gets its data from a XMLInputStream object that is
	 * created (internally) to work with this source. 
	 *
	 * getObject() is implemented with EXPAT if HAVE_EXPAT is 1, or with
	 * Xerces if HAVE_XERCES is 1.
	 */
	class XMLInputSource
		{
		public:
			XMLInputSource(){;}

			/** virtual dtor */
			virtual ~XMLInputSource(){;}
			
			/**  Virtual input source interface */
			virtual TSFSmartPtr<XMLInputStream> stream() const = 0 ;
			
			/** get an object by invoking the TreeBuildingXMLHandler on the
			 * input data */
			XMLObject getObject() const ;

		};

}
#endif

