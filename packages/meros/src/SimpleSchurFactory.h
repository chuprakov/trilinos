#ifndef SIMPLESCHURFACTORY_H
#define SIMPLESCHURFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"
#include "TSFOperatorSourceBase.h"
#include "RightBlockNSOperatorSource.h"
#include "SimpleOperatorSource.h"
#include "SchurFactoryBase.h"
#include "SchurFactory.h"
#include "SimpleOperatorSource.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{

  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * SimpleSchurFactory
   *
   */

  class SimpleSchurFactory : public SchurFactoryBase
    {
    public:
      /** ctor */
      SimpleSchurFactory(TSFLinearSolver& schurSolver);

      /** virtual destructor */
      virtual ~SimpleSchurFactory(){;}

      /** get Schur complement inverse approximation */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& diagfac) const;

      /** write to a string */
      virtual string toString() const ;

    private:
      TSFLinearSolver schurSolver_;
    };

  /** \relates DiagSchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const SimpleSchurFactory& x);

  /** \relates DiagSchurFactory
   * write to a string
   */
  string toString(const SimpleSchurFactory& x);
}


#endif
