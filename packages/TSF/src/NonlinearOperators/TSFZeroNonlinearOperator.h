#ifndef TSFZERONONLINEAROPERATOR_H
#define TSFZERONONLINEAROPERATOR_H

#include "TSFDefs.h"
#include "TSFNonlinearOperatorBase.h"

namespace TSF
{


  /** \ingroup NonlinearOperatorSubtypes
   * TSFZeroNonlinearOperator is a nonlinear operator that returns
   * the zero vector in the range space.
   */

  class TSFZeroNonlinearOperator : public TSFNonlinearOperatorBase
    {
    public:
      /** construct with domain and range spaces  */
      TSFZeroNonlinearOperator(const TSFVectorSpace& domain,
                               const TSFVectorSpace& range);

      /** the usual virtual dtor */
      virtual ~TSFZeroNonlinearOperator(){;}

      /** apply() returns a zero vector for all input */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /** the derivative is the zero operator */
      virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

      /** write to a stream */
      virtual void print(ostream& os) const ;
    protected:
    };
}

#endif
