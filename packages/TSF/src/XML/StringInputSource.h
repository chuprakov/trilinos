#ifndef STRINGINPUTSOURCE_H
#define STRINGINPUTSOURCE_H

#include "TSFDefs.h"

#include "XMLInputSource.h"
#include <string>

namespace TSF
{

  using std::string;

  /** \ingroup XML
   * StringInputSource reads XML from a String
   */

  class StringInputSource : public XMLInputSource
    {
    public:
      /** ctor */
      StringInputSource(const string& text);
      /** virtual dtor */
      virtual ~StringInputSource(){;}

      /** create a StringInputStream */
      virtual TSFSmartPtr<XMLInputStream> stream() const;

    private:
      string text_;
    };


}
#endif

