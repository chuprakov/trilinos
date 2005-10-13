#ifndef CGSOLVERVEH_H
#define CGSOLVERVEH_H

#include "TSFConfig.h"
#include "TSFLinearSolverBase.h"
#include "TSFPreconditionerFactory.h"
#include "TSFHashtable.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Representation-independent CG solver. This is just for testing,
   * so there is no provision for preconditioning.
   */

  class CGSolverVEH : public TSFLinearSolverBase
    {
    public:
      /** construct a CG solver with
       * parameters for max iterations and convergence tolerance */
      CGSolverVEH(const TSFReal& tol = 1.0e-08, 
		  int maxIters = 300);

      /** construct a preconditioned CG solver with
       * parameters for max iterations and convergence tolerance */
      CGSolverVEH(const TSFPreconditionerFactory& pf,
                  const TSFReal& tol = 1.0e-08, 
		  int maxIters = 300);

      /** TUVD */
      virtual ~CGSolverVEH();


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

/*       /\* Initial solve call -- set up x0, etc. *\/ */
/*       virtual bool solveStart(const TSFLinearOperator& op, */
/*                               const TSFVector& rhs, */
/*                               TSFVector& soln) const ; */

/*      virtual bool solveRestart(TSFVector& soln) const ; */

      virtual bool solveUnpreconditioned(const TSFLinearOperator& op,
					 const TSFVector& rhs,
					 TSFVector& soln) const ;


      TSFReal tol_;
      int maxIters_;
      /* factory to build a preconditioner */
      TSFPreconditionerFactory preconditionerFactory_;
    };

}


#endif
