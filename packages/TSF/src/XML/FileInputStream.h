#ifndef FILEINPUTSTREAM_H
#define FILEINPUTSTREAM_H

#include "TSFDefs.h"
#include "XMLInputStream.h"
#include <string>
#include <cstdio>

namespace TSF
{
  using std::string;

  /** \ingroup XML
   * Input stream for reading an entire document from a file. This
   * is a low-level object and should not be needed at the user level.
   * FileInputSource is the user-level object.
   */

  class FileInputStream : public XMLInputStream
    {
    public:
      /** construct with a filename */
      FileInputStream(const string& filename);
      /** TUVD */
      virtual ~FileInputStream() {;}

      /** read up to maxToRead bytes */
      virtual unsigned int readBytes(unsigned char* const toFill,
                                     const unsigned int maxToRead);

    private:
      FILE* file_;
    };
}
#endif

