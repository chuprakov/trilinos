#ifndef TREEBUILDINGXMLHANDLER_H
#define TREEBUILDINGXMLHANDLER_H

#include "TSFDefs.h"
#include "XMLObject.h"
#include "TSFStack.h"


namespace TSF
{

  /** \ingroup XML
   * Use TreeBuildingXMLHandler assembles
   * a XMLObject from your XML input.
   */

  class TreeBuildingXMLHandler
    {
    public:
      /** empty ctor only */
      TreeBuildingXMLHandler();

      /** retrieve the entire XML tree */
      const XMLObject& getObject() const {return root_.getChild(0);}

      /** Process character data */
      void characters(const string& chars, const unsigned int length);

      /** Receive notification of the end of an element */
      void endElement(const string& tag);

      /** Receive notification of the start of an element */
      void startElement(const string& tag,
                        const TSFHashtable<string,string>& attributes);

    private:
      XMLObject root_;
      XMLObject current_;
      TSFStack<XMLObject> path_;
    };
}


#endif

