#ifndef BJBRIGHTOPERATORSOURCE_H
#define BJBRIGHTOPERATORSOURCE_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * BJBRightOperatorSource
   * X = C - B D^{-1} B^T, where D = diag(conv-diff op)
   *
   */

  // Need a better name for this class.

  class BJBRightOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      // BJBRightOperatorSource(TSFLinearOperator& S);

      /** ctor (we make Ap = C from saddle matrix S) */
      BJBRightOperatorSource(TSFLinearOperator& S);

      /** ctor (we are given Ap) */
      BJBRightOperatorSource(TSFLinearOperator& S,
                             TSFLinearOperator& Ap);


      /** virtual destructor */
      virtual ~BJBRightOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getAp() const;


      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      mutable TSFLinearOperator Ap_;
      mutable bool hasAp_;
    };


}


#endif
