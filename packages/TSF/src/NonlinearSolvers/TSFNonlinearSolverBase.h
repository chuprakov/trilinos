#ifndef TSFNONLINEARSOLVERBASE_H
#define TSFNONLINEARSOLVERBASE_H

#include "TSFDefs.h"
#include "TSFLinearSolver.h"
#include "TSFNonlinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup NonlinearSolverSubtypes
   *
   */

  class TSFNonlinearSolverBase
    {
    public:
      /** construct with a maximum number of iterations */
      TSFNonlinearSolverBase(int maxIters);

      /** TUVD */
      virtual ~TSFNonlinearSolverBase(){;}

      /**
       * Solve the system F(x)==0 for x, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFNonlinearOperator& op,
                         TSFVector& soln) const = 0 ;

    protected:
      int maxIters_;
    };

}


#endif
