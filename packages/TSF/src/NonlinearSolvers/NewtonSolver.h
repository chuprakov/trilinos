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
      /** */
      NewtonSolver(const TSFLinearSolver& linearSolver,
                   int maxIters,
                   const double& stepTol,
                   const double& funcTol,
                   const double& low, const double& high);

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
      bool hasBounds() const {return hasBounds_;}

      double distanceToBoundary(const TSFVector& x0, const TSFVector& delta) const ;

      TSFLinearSolver linearSolver_;
      double stepTol_;
      double funcTol_;
      bool hasBounds_;
      double low_;
      double high_;
    };

}


#endif
