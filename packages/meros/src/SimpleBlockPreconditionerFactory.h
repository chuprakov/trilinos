#ifndef SIMPLEBLOCKPRECONDITIONERFACTORY_H
#define SIMPLEBLOCKPRECONDITIONERFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFLinearSolver.h"
#include "TSFArray.h"
#include "SchurFactoryBase.h"
#include "SchurFactory.h"
#include "SimpleOperatorSource.h"
#include <iostream>
#include <string>

namespace SPP
{
  using namespace TSF;
  /** \ingroup Preconditioner
   *  factory for block preconditioners for Simple
   */
  class SimpleBlockPreconditionerFactory : public TSFPreconditionerFactoryBase
    {
    public:
      /** constructor for preconditioner using right LDU factors only */
      SimpleBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac);

      SimpleBlockPreconditionerFactory(const TSFLinearSolver& Fsolver, const SchurFactory& sfac, const TSFLinearSolver& SchurSolver);

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFLinearOperator& op) const;

      /** create a concrete preconditioner */
      virtual TSFPreconditioner createPreconditioner(const TSFOperatorSource& saddleOpSrc) const;

      /** write to a string */
      virtual string toString() const;

    private:
      TSFLinearSolver Fsolver_;
      TSFLinearSolver Schursolver_;
      SchurFactory sfac_;
      TSFLinearOperator op_;
    };

}
#endif
