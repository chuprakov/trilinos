#ifndef BJBSCHURFACTORY_H
#define BJBSCHURFACTORY_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"
#include "TSFLinearSolver.h"
#include "TSFOperatorSourceBase.h"
#include "RightBlockNSOperatorSource.h"
#include "BJBRightOperatorSource.h"
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
   * BJBSchurFactory
   *
   */

  class BJBSchurFactory : public SchurFactoryBase
    {
    public:
      /** ctor */
      BJBSchurFactory(TSFLinearSolver& ApSolver);

      //      /** ctor (with Mp solver) */
      //      BJBSchurFactory(TSFLinearSolver& ApSolver,
      //                            TSFLinearSolver& MpSolver);

      // maybe have a default that takes no solvers?

      /** virtual destructor */
      virtual ~BJBSchurFactory(){;}

      /** get Schur complement inverse approximation */
      virtual TSFLinearOperator
        getSchurInvApprox(const TSFOperatorSource& bjbfac) const;

      /** write to a string */
      virtual string toString() const ;

    private:
      TSFLinearSolver ApSolver_;
    };

  /** \relates BJBSchurFactory
   * write to an output stream
   */
  ostream& operator<<(ostream& os, const BJBSchurFactory& x);

  /** \relates BJBSchurFactory
   * write to a string
   */
  string toString(const BJBSchurFactory& x);



}


#endif
