#ifndef KAYLOGHINSCHURFACTORY_H
#define KAYLOGHINSCHURFACTORY_H

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
   * KayLoghinSchurFactory
   *
   */

  class KayLoghinSchurFactory : public SchurFactoryBase
    {
    public:
      /** ctor */
      KayLoghinSchurFactory(TSFLinearSolver& ApSolver);

      //      /** ctor (with Mp solver) */
      //      KayLoghinSchurFactory(TSFLinearSolver& ApSolver,
      //                            TSFLinearSolver& MpSolver);

      // maybe have a default that takes no solvers?

      /** virtual destructor */
      virtual ~KayLoghinSchurFactory(){;}

      /** get Schur complement inverse approximation */
      virtual TSFLinearOperator getSchurInvApprox(const TSFOperatorSource& klfac) const;

      /** write to a string */
      virtual string toString() const ;

    private:
      TSFLinearSolver ApSolver_;
    };

  /** \relates KayLoghinSchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const KayLoghinSchurFactory& x);

  /** \relates KayLoghinSchurFactory
   * write to a string
   */
  string toString(const KayLoghinSchurFactory& x);



}


#endif
