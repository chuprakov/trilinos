/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFACCESSIBLEVECTOR_HPP
#define TSFACCESSIBLEVECTOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreTypes.hpp"



#ifndef DOXYGEN_DEVELOPER_ONLY

namespace TSFExtended
{
  using TSFCore::Index;

  /**
   * TSFExtended::AccessibleVector defines an interface through which
   * elements for a vector can be accessed. Element loading is occasionally
   * used by application codes in probing results vectors, 
   * but should rarely be used by high-performance solver codes; this 
   * capability is therefore in TSFExtended rather than TSFCore.
   *
   * @author Kevin Long (krlong@sandia.gov)
   */
  template <class Scalar>
  class AccessibleVector 
    {
    public:
      /** virtual dtor */
      virtual ~AccessibleVector() {;}

      /** get the element at the given global index */
      virtual const Scalar& getElement(Index globalIndex) const = 0 ;
    };
}

#endif  /* DOXYGEN_DEVELOPER_ONLY */

#endif
