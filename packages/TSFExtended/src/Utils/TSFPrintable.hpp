/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFPRINTABLE_HPP
#define TSFPRINTABLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"

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

#endif
