#ifndef TSFCOMPOSEDOPERATOR_H
#define TSFCOMPOSEDOPERATOR_H

#include "TSFDefs.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"

namespace TSF
{


  /** \ingroup LinearOperatorSubtypes
   * TSFComposedOperator is a composition of two linear operators.
   */

  class TSFComposedOperator : public TSFLinearOperatorBase
    {
    public:
      /** construct with a pair of operators and a boolean to indicate
       * if addition or substraction is to be performed. */
      TSFComposedOperator(const TSFLinearOperator& left,
                          const TSFLinearOperator& right);

      /* the usual virtual dtor */
      virtual ~TSFComposedOperator(){;}

      /** apply operator to a vector in the domain space, returning a vector
       * in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;
      /** apply the operator's inverse to a vector. If the composed operator
       * is A*B, the inverse will be inv(B)*inv(A). */
      void applyInverse(const TSFVector& in, TSFVector& out) const ;

      /** apply the operator's adjoint to a vector. If the composed operator
       * is A*B, the adjoint will be adj(B)*adj(A). */
      void applyAdjoint(const TSFVector& in, TSFVector& out) const ;

    protected:
      TSFLinearOperator left_;
      TSFLinearOperator right_;
    };
}

#endif
