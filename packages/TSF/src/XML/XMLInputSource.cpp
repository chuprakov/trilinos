#include "XMLInputSource.h"
#include "TreeBuildingXMLHandler.h"


#if HAVE_XERCES
#include "XercesHandlerAdapter.h"
#include "XercesInputSourceAdapter.h"
#include "XercesInputStreamAdapter.h"
#include <util/PlatformUtils.hpp>
#elif HAVE_EXPAT
#include "ExpatHandlerAdapter.h"
#define EXPAT_BUFSIZE 8192
#endif


using namespace TSF;

XMLObject XMLInputSource::getObject() const
{
#if HAVE_XERCES

	static bool first = true;
	if (first)
		{
			XMLPlatformUtils::Initialize();
			first = false;
		}

	SAXParser parser;
	XercesHandlerAdapter handler(new TreeBuildingXMLHandler());
	XercesInputSourceAdapter inputSource(this);

	parser.setDocumentHandler(&handler);

	try 
		{
			parser.parse(inputSource);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in SAX parsing");
		}
	return handler.getObject();

#elif HAVE_EXPAT

	TSFSmartPtr<TreeBuildingXMLHandler> handler = new TreeBuildingXMLHandler();

	XML_Parser parser = XML_ParserCreate(NULL);

	XML_SetElementHandler(parser, expatStartElementHandler, 
												expatEndElementHandler);

	XML_SetCharacterDataHandler(parser, expatCharacterDataHandler);

	XML_SetUserData(parser, (void*) &(*handler));

	TSFSmartPtr<XMLInputStream> s = stream();

	bool done = false;
	unsigned int bufsize = EXPAT_BUFSIZE;
	unsigned char buf[EXPAT_BUFSIZE];

	while (!done)
		{
			unsigned int nRead = s->readBytes(buf, bufsize);
			if (nRead < bufsize) 
				{
					done = true;
				}
			XML_Parse(parser, (char*) buf, bufsize, done);
		}

	return handler->getObject();

#else
	TSFError::raise("XMLInputSource::getObject() - no XML parser installed");
	return XMLObject(); // -Wall
#endif

}
			
			
