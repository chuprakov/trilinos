#ifndef DIAGSCHURFACTORY_H
#define DIAGSCHURFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"
#include "TSFOperatorSourceBase.h"
#include "RightBlockNSOperatorSource.h"
#include "KayLoghinRightOperatorSource.h"
#include "SchurFactoryBase.h"
#include "SchurFactory.h"
#include "TSFArray.h"
#include <iostream>
#include <string>

namespace SPP
{

  using namespace TSF;
  using std::string;
  using std::ostream;


  /** \ingroup Preconditioner
   * DiagSchurFactory
   *
   */

  class DiagSchurFactory : public SchurFactoryBase
    {
    public:
      /** ctor */
      DiagSchurFactory(TSFLinearSolver schurSolver);

      /** virtual destructor */
      virtual ~DiagSchurFactory(){;}

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
  ostream& operator<<(ostream& os, const DiagSchurFactory& x);

  /** \relates DiagSchurFactory
   * write to a string
   */
  string toString(const DiagSchurFactory& x);



}


#endif
