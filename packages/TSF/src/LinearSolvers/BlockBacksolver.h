#ifndef BLOCKBACKSOLVER_H
#define BLOCKBACKSOLVER_H

#include "TSFDefs.h"
#include "TSFArray.h"
#include "TSFLinearSolver.h"

namespace TSF
{
  /** block solver for NxN block upper triangular systems. */

  class BlockBacksolver : public TSFLinearSolverBase
    {
    public:
      /** */
      BlockBacksolver(const TSFLinearSolver& s);

      /** */
      BlockBacksolver(const TSFLinearSolver& s0, const TSFLinearSolver& s1);

      /** */
      BlockBacksolver(const TSFLinearSolver& s0,
                      const TSFLinearSolver& s1,
                      const TSFLinearSolver& s2);

      /** */
      BlockBacksolver(const TSFArray<TSFLinearSolver>& s);

      /** the usual virtual dtor */
      virtual ~BlockBacksolver() {;}

      /** you guessed it, solve the system  */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;
    private:
      TSFArray<TSFLinearSolver> solvers_;
    };

}

#endif
