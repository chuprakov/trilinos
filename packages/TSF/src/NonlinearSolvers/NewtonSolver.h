#ifndef NEWTONSOLVER_H
#define NEWTONSOLVER_H

#include "TSFDefs.h"
#include "TSFNonlinearSolver.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup NonlinearSolverSubtypes
   *
   */

  class NewtonSolver : public TSFNonlinearSolverBase
    {
    public:
      /** */
      NewtonSolver(const TSFLinearSolver& linearSolver,
                   int maxIters,
                   const double& stepTol,
                   const double& funcTol);

      /** TUVD */
      virtual ~NewtonSolver();

      /**
       * Solve the system F(x)==0 for x, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFNonlinearOperator& op,
                         TSFVector& soln) const ;

      /**
       * Solve the system F(x)==0 for x, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFNonlinearOperatorBase& op,
                         TSFVector& soln) const ;

    private:
      TSFLinearSolver linearSolver_;
      double stepTol_;
      double funcTol_;
    };

}


#endif
