/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFDESCRIBABLE_HPP
#define TSFDESCRIBABLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtended
{
  /**
   * Describable defines an interface for writing short descriptive strings
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  class Describable
    {
    public:
      /** Return a brief descriptive string */
      virtual string describe() const = 0 ;
    };
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
