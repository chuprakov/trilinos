#ifndef SCHURCOMPLEMENTSOLVER_H
#define SCHURCOMPLEMENTSOLVER_H

#include "TSFDefs.h"
#include "TSFLinearSolverBase.h"
#include "TSFLinearSolver.h"
#include "TSFPreconditionerFactory.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Representation-independent Schur complement solver. This solves a block
   * linear system of the form [A, B ; C, D]*[x, y]' = [a, b]'
   * by reducing to the triangular system
   *    [ A, B ; 0, D - C A^(-1)B ] * [x, y]'   [a,  b - C A^(-1) a ]'
   *
   * and backsolving. Two solvers are needed: an "inner" solver to carry out the
   * solves on A needed in the application of A^(-1), and an "outer" solver
   * to carry out the solve on the schur complement D - C A^(-1) B.
   */

  class SchurComplementSolver : public TSFLinearSolverBase
    {
    public:
      SchurComplementSolver(const TSFLinearSolver& innerSolver,
                            const TSFLinearSolver& outerSolver);

      virtual ~SchurComplementSolver();

      /**
       * Solve the system with the given RHS, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed. The input operator is expected to be a 2x2 block operator,
       * and the input rhs is expected to be a
       */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;
    private:

      TSFLinearSolver innerSolver_;

      TSFLinearSolver outerSolver_;
    };

}


#endif
