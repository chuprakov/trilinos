#ifndef CGSOLVER_H
#define CGSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Representation-independent CG solver. This is just for testing,
   * so there is no provision for preconditioning.
   */

  class CGSolver : public TSFLinearSolverBase
    {
    public:
      /** construct a CG solver to work with a given operator, with
       * parameters for max iterations and convergence tolerance */
      CGSolver(const TSFParameterList& params);

      /**
       * Solve the system with the given RHS, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;

      /** */
      virtual TSFParameterList defaultParameters() const ;
    private:
      TSFReal tol_;
      int maxIters_;
    };

}


#endif
