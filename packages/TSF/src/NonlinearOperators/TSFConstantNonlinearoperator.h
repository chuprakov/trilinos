#ifndef TSFCONSTANTNONLINEAROPERATOR_H
#define TSFCONSTANTNONLINEAROPERATOR_H

#include "TSFDefs.h"
#include "TSFNonlinearOperatorBase.h"

namespace TSF
{


  /** \ingroup NonlinearOperatorSubtypes
   * TSFConstantNonlinearOperator is a nonlinear operator that returns
   * a constant vector in the range space.
   */

  class TSFConstantNonlinearOperator : public TSFLinearOperatorBase
    {
    public:
      /** construct with a vector giving the constant return value */
      TSFConstantNonlinearOperator(const TSFVectorSpace& domain,
                                   const TSFVector& value);

      /** the usual virtual dtor */
      virtual ~TSFConstantNonlinearOperator(){;}

      /** apply operator to a vector in the domain space, returning a vector
       * in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /** get a linear operator representing the derivative */
      virtual TSFLinearOperator derivative() const ;

    protected:
      TSFVector value_;
    };
}

#endif
