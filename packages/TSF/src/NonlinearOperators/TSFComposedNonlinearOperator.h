#ifndef TSFCOMPOSEDNONLINEAROPERATOR_H
#define TSFCOMPOSEDNONLINEAROPERATOR_H

#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFNonlinearOperator.h"

namespace TSF
{


  /** \ingroup NonlinearOperatorSubtypes
   * TSFComposedNonlinearOperator is the composition of two nonlinear
   * operators, left() and right(), giving y = left(right(x)).
   */

  class TSFComposedNonlinearOperator : public TSFNonlinearOperatorBase
    {
    public:
      /** construct with a pair of operators.  */
      TSFComposedNonlinearOperator(const TSFNonlinearOperator& left,
                                   const TSFNonlinearOperator& right);
      /** the usual virtual dtor */
      virtual ~TSFComposedNonlinearOperator(){;}

      /** apply operator to a vector in the domain space, returning a vector
       * in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /** get a linear operator representing the derivative */
      virtual TSFLinearOperator derivative(const TSFVector& evalPt) const ;

      /** write to a stream */
      virtual void print(ostream& os) const ;
    protected:
      TSFNonlinearOperator left_;
      TSFNonlinearOperator right_;
    };
}

#endif
