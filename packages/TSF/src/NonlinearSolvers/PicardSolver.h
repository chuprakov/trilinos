#ifndef PICARDSOLVER_H
#define PICARDSOLVER_H

#include "TSFDefs.h"
#include "TSFNonlinearSolver.h"
#include "TSFConvergenceTest.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup NonlinearSolverSubtypes
   *
   */

  class PicardSolver : public TSFNonlinearSolverBase
    {
    public:
      /** */
      PicardSolver(int maxIters,
                   const double& tol,
                   const double& relax = 0.0);

      /** TUVD */
      virtual ~PicardSolver();

      /**
       * Solve the system F(x)==x for x,
       */
      virtual bool solve(const TSFNonlinearOperator& op,
                         TSFVector& soln) const ;
      /**
       * Solve the system F(x)==x for x,
       */
      virtual bool solve(const TSFNonlinearOperatorBase& op,
                         TSFVector& soln) const ;
    private:
      double tol_;
      double relax_;
    };

}


#endif
