/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFVECTOR_HPP
#define TSFVECTOR_HPP

#include "TSFConfigDefs.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   *
   */
  template <class Scalar>
  class VectorSpace : public Handle<TSFCore::VectorSpace<Scalar> >
  {
  public:
    /** */
    VectorSpace(TSFCore::VectorSpace* ptr);
    
    /** */
    Vector<Scalar> createMember() const ;
  };
  

}


#endif
