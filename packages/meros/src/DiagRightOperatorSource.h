#ifndef DIAGRIGHTOPERATORSOURCE_H
#define DIAGRIGHTOPERATORSOURCE_H

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
   * DiagRightOperatorSource
   * X = C - B D^{-1} B^T, where D = diag(conv-diff op)
   *
   */

  // Need a better name for this class.

  class DiagRightOperatorSource : public RightBlockNSOperatorSource
    {
    public:
      /** ctor */
      // DiagRightOperatorSource(TSFLinearOperator& S);

      /** ctor (we make D_inv = inv(diag(F)) from F) */
      DiagRightOperatorSource(TSFLinearOperator& S);

      /** ctor (we are given D = diag(F) */
      DiagRightOperatorSource(TSFLinearOperator& S,
                              TSFLinearOperator& Dinv);


      /** virtual destructor */
      virtual ~DiagRightOperatorSource(){;}

      /** get a concrete linear operator */
      TSFLinearOperator getOp() const;

      /** get pressure Poisson operator */
      TSFLinearOperator getDinv() const;


      /** write to a string */
      string toString() const ;

    private:
      TSFLinearOperator S_;
      mutable TSFLinearOperator Dinv_;
      mutable bool hasDinv_;
    };


}


#endif
