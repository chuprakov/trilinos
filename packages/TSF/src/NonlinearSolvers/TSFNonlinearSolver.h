#ifndef TSFNONLINEARSOLVER_H
#define TSFNONLINEARSOLVER_H

#include "TSFDefs.h"
#include "TSFNonlinearSolverBase.h"
#include "TSFLinearSolver.h"
#include "TSFSmartPtr.h"

namespace TSF
{


  /** \ingroup Solvers
   *
   */

  class TSFNonlinearSolver
    {
    public:
      /** empty ctor */
      TSFNonlinearSolver();
      /** construct with a ptr to the base class */
      TSFNonlinearSolver(TSFNonlinearSolverBase* ptr);

      /** \name Solve methods */
      //@{
      /**
       * Solve the system F(x)==0 for x, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      bool solve(const TSFNonlinearOperator& op,
                 TSFVector& soln) const ;

      /**
       * Solve the system F(x)==y for x, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      bool solve(const TSFNonlinearOperator& op,
                 const TSFVector& y,
                 TSFVector& soln) const ;

      //@}
      /** */
      bool isNull() const ;

    private:
      TSFSmartPtr<TSFNonlinearSolverBase> ptr_;
    };

}


#endif
