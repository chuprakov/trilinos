#ifndef GMRESSOLVER_H
#define GMRESSOLVER_H

#include "TSFDefs.h"
#include "TSFLinearSolverBase.h"
#include "TSFPreconditionerFactory.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Representation-independent GMRES solver.
   */

  class GMRESSolver : public TSFLinearSolverBase
    {
    public:
      /** construct a GMRES solver to work with a given operator, with
       * parameters for max iterations, convergence tolerance, and restart value m */
      GMRESSolver(const TSFReal& tol = 1.0e-08, int maxIters = 300, int m = 50);

      /** construct a GMRES solver to work with a given operator, with
       * parameters for max iterations, convergence tolerance, and restart value m */
      GMRESSolver(const TSFPreconditionerFactory& pf,
                  const TSFReal& tol = 1.0e-08, int maxIters = 300, int m = 50);

      /** TUVD */
      virtual ~GMRESSolver(){;}
      /** what is that??? -- veh */


      /**
       * Solve the system with the given RHS, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;
    private:
      virtual bool solveUnpreconditioned(const TSFLinearOperator& op,
                                         const TSFVector& rhs,
                                         TSFVector& soln) const ;


      /* residual convergence tolerance */
      TSFReal tol_;
      /* maximum total number of iterations */
      int maxIters_;
      /* restart value */
      int m_;
      /* factory to build a preconditioner */
      TSFPreconditionerFactory preconditionerFactory_;
    };

}


#endif
