/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFEXPLICITLYTRANSPOSEABLEOP_HPP
#define TSFEXPLICITLYTRANSPOSEABLEOP_HPP

#include "TSFConfigDefs.hpp"


namespace TSFExtended
{

  /** 
   * Base interface for operators whose transpose can be formed
   * explicitly. Note: if a transpose is formed, it will be an independent
   * object in the sense that any changes made to either the original object
   * or the transpose are not propagated to the other.
   */
  template <class Scalar>
  class ExplicitlyTransposeableOp 
    {
    public:
      /** Virtual dtor */
      virtual ~ExplicitlyTransposeableOp(){;}

      /** Form a new object that is the transpose of this operator */
      virtual LinearOperator<Scalar> formTranspose() const = 0 ;
      
    private:
    };
}


#endif
