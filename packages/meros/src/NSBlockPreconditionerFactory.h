#ifndef NSBLOCKPRECONDITIONERFACTORY_H
#define NSBLOCKPRECONDITIONERFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFLinearSolver.h"
#include "TSFArray.h"
#include "SchurFactoryBase.h"
#include "SchurFactory.h"
#include "KayLoghinRightOperatorSource.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  /** \ingroup Preconditioner
   *  factory for block preconditioners for NS
   */
  class NSBlockPreconditionerFactory : public TSFPreconditionerFactoryBase
    {
    public:
      /** constructor for preconditioner using right LDU factors only */
      NSBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac);

      // /** constructor for preconditioner using left and right LDU factors*/
      // NSBlockPreconditionerFactory(const MomentumFactory& mfac,
      //                             const SchurFactory& sfac,
      //                             const ProjectionFactory& projfac);

      //      /** virtual destructor */
      //      virtual ~NSBlockPreconditionerFactory(){;}

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& op) const;

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFOperatorSource& saddleOpSrc) const;

      /** write to a string */
      virtual string toString() const;

    private:
      TSFLinearSolver Fsolver_;
      SchurFactory sfac_;
      TSFLinearOperator op_;
    };

}
#endif
