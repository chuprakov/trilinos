#ifndef TSFZEROOPERATOR_H
#define TSFZEROOPERATOR_H

#include "TSFDefs.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFArray.h"

namespace TSF
{




  /** \ingroup LinearOperatorSubtypes
   * TSFZeroOperator is the zero operator, which maps any vector in the
   * domain space to the zero vector in the range space. The adjoint of the
   * zero operator maps any vector in the range space (the domain of the
   * adjoint) to the zero vector in the domain space (the range of the
   * adjoint). The inverse of the zero operator is undefined, so calling
   * applyInverse() results in an error.
   */

  class TSFZeroOperator : public TSFLinearOperatorBase
    {
    public:
      /** Construct with the domain and range spaces */
      TSFZeroOperator(const TSFVectorSpace& domain,
                      const TSFVectorSpace& range);

      /** the usual virtual dtor */
      virtual ~TSFZeroOperator(){;}

      /** identify self as a zero operator */
      virtual bool isZeroOperator() const {return true;}

      /** apply returns a zero vector in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /** applyAdjoint returns a zero vector in the domain space */
      virtual void applyAdjoint(const TSFVector& in,
                                TSFVector& out) const ;

      /** calling applyInverse is an error */
      virtual void applyInverse(const TSFVector& in,
                                TSFVector& out) const ;

      /**  get row  */
      virtual void getRow(int row, TSFArray<int>& indices,
                          TSFArray<TSFReal>& values) const;


      /**
       * Write to a stream
       */
      virtual void print(ostream& os) const ;
    protected:
    };
}

#endif
