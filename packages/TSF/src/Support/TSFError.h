#ifndef TSFERROR_H
#define TSFERROR_H

#include "TSFDefs.h"
#include <string>
#include <stdexcept>
#include "TSFRaiseHandlerBase.h"
#include "TSFSmartPtr.h"

namespace TSF
{
  using std::exception;
  using std::string;

  /** \ingroup ErrorHandling
   *
   */
  class TSFError
    {
    public:

      /** raise() is to be called when an error is detected. */
      static void raise(const string& msg);

      /** */
      static void trace(const exception& what, const string& where);

      /** setRaiseHandler() is a hook for installation of a custom raise
       * handler*/
      static void setRaiseHandler(TSFRaiseHandlerBase* handler);

    private:
      static TSFSmartPtr<TSFRaiseHandlerBase> handler_;
    };
}
#endif
