#ifndef TSFRAISEHANDLERBASE_H
#define TSFRAISEHANDLERBASE_H

#include "TSFDefs.h"
#include <string>

namespace TSF
{
  using std::string;

  /** \ingroup ErrorHandling
   * TSFRaiseHandlerBase is a base class for objects that specify what is to
   * happen upon a call to TSFError::raise(). When TSFError::raise() is called,
   * the handleRaise() method of a raise handler subclass will be called.
   * The behavior can be customized by supplying a custom raise hander.
   * By default, the TSFDefaultRaiseHandler is used and the handleRaise()
   * method throws an exception.
   */
  class TSFRaiseHandlerBase
    {
    public:
      /** empty ctor */
      TSFRaiseHandlerBase(){;}

      /** TUVD */
      virtual ~TSFRaiseHandlerBase(){;}

      /** the handleRaise() method is called inside TSFError::raise() */
      virtual void handleRaise(const char* msg) = 0 ;

    private:

    };
}
#endif
