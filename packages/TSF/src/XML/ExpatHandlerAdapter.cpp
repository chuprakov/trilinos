#include "ExpatHandlerAdapter.h"

#if HAVE_EXPAT

#include "TreeBuildingXMLHandler.h"

using namespace TSF;

void expatStartElementHandler(void* handler, 
															const XML_Char* name, 
															const XML_Char** attr)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	string tag = name;
	TSFHashtable<string, string> attributes;
	
	/* the attribute data is stored in a C array of C strings, in order 
	 * {key1, val1, key2, val2, ...}. */

	for (int i=0; attr[i] != 0; i+=2)
		{
			string key = attr[i];
			string val = attr[i+1];
			attributes.put(key, val);
		}

	h->startElement(tag, attributes);
}

void expatEndElementHandler(void* handler, 
														const XML_Char* name)
{
	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	
	string tag = name;

	h->endElement(tag);
}

void expatCharacterDataHandler(void* handler, 
															 const XML_Char* s,
															 int len)
{
	char* str = new char[len+1];
	strncpy(str, s, len);

	string chars = str;

	TreeBuildingXMLHandler* h = (TreeBuildingXMLHandler*) handler;
	h->characters(chars, len);
}

#endif
