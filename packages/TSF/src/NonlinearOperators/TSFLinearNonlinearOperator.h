#ifndef TSFLINEARNONLINEAROPERATOR_H
#define TSFLINEARNONLINEAROPERATOR_H

#include "TSFDefs.h"
#include "TSFNonlinearOperatorBase.h"
#include "TSFLinearOperatorBase.h"

namespace TSF
{


  /** \ingroup NonlinearOperatorSubtypes
   * TSFLinearNonlinearOperator is a linear operator represented as
   * a nonlinear operator object.
   */

  class TSFLinearNonlinearOperator : public TSFNonlinearOperatorBase
    {
    public:
      /** construct with a vector giving the constant return value */
      TSFLinearNonlinearOperator(const TSFLinearOperator& op)
        : TSFNonlinearOperatorBase(op.domain(), op.range()), op_(op) {;}

      /** the usual virtual dtor */
      virtual ~TSFLinearNonlinearOperator(){;}

      /** apply operator to a vector in the domain space, returning a vector
       * in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const {op_.apply(in, out);}

      /** get a linear operator representing the derivative */
      virtual TSFLinearOperator derivative(const TSFVector& /* evalPt */) const
        {return op_;}

      /** write to a stream */
      virtual void print(ostream& os) const {os << op_;}
    protected:
      TSFLinearOperator op_;
    };
}

#endif
