#ifndef FILEINPUTSOURCE_H
#define FILEINPUTSOURCE_H

#include "TSFDefs.h"

#include "XMLInputSource.h"
#include <string>

namespace TSF
{

  /** \ingroup XML
   * Input source that reads XML from a file.
   */

  class FileInputSource : public XMLInputSource
    {
    public:
      /** ctor */
      FileInputSource(const string& filename);

      /** virtual dtor */
      virtual ~FileInputSource(){;}

      /** create a FileInputStream */
      virtual TSFSmartPtr<XMLInputStream> stream() const;

    private:
      string filename_;
    };

}
#endif

