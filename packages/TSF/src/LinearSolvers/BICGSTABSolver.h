#ifndef BICGSTABSOLVER_H
#define BICGSTABSOLVER_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"
#include "TSFPreconditionerFactory.h"
#include "TSFTimeMonitor.h"
#include "TSFParameterListImplem.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Representation-independent BICGSTAB solver.
   */

  class BICGSTABSolver : public TSFLinearSolverBase
    {
    public:
      /** construct a BICGSTAB solver with
       * parameters for max iterations and convergence tolerance */
      BICGSTABSolver(const TSFParameterList& params = TSFParameterList());

      /** construct a preconditioned BICGSTAB solver, with
       * parameters for max iterations and convergence tolerance */
      BICGSTABSolver(const TSFPreconditionerFactory& pf,
                     const TSFParameterList& params = TSFParameterList());

      /** Deprecated: contruct with explicit arguments for tolerance and maxiters.
       * New code should use parameters. */
      BICGSTABSolver(const double& tol, int maxiters);

      /** Deprecated: contruct with explicit arguments for tolerance and maxiters.
       * New code should use parameters. */
      BICGSTABSolver(const TSFPreconditionerFactory& pf,
                     const double& tol, int maxiters);

      /** TUVD */
      virtual ~BICGSTABSolver();

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

      virtual bool solveUnpreconditioned(const TSFLinearOperator& op,
                                         const TSFVector& rhs,
                                         TSFVector& soln) const ;
      /* factory to build a preconditioner */
      TSFPreconditionerFactory preconditionerFactory_;
    };

}


#endif
