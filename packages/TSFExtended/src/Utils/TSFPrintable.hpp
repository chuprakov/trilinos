/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFPRINTABLE_HPP
#define TSFPRINTABLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"


#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtended
{
  /**
   * Printable defines an interface for writing an object to a stream. 
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  class Printable
    {
    public:
      virtual void print(ostream& os) const = 0 ;
    };
}


#endif  /* DOXYGEN_DEVELOPER_ONLY */
#endif
