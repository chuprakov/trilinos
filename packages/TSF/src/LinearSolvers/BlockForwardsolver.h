#ifndef BLOCKFORWARDSOLVER_H
#define BLOCKFORWARDSOLVER_H

#include "TSFDefs.h"
#include "TSFArray.h"
#include "TSFLinearSolver.h"

namespace TSF
{
  /** block solver for NxN block upper triangular systems. */

  class BlockForwardsolver : public TSFLinearSolverBase
    {
    public:
      /** */
      BlockForwardsolver(const TSFLinearSolver& s);

      /** */
      BlockForwardsolver(const TSFLinearSolver& s0, const TSFLinearSolver& s1);

      /** */
      BlockForwardsolver(const TSFLinearSolver& s0,
                      const TSFLinearSolver& s1,
                      const TSFLinearSolver& s2);

      /** */
      BlockForwardsolver(const TSFArray<TSFLinearSolver>& s);

      /** the usual virtual dtor */
      virtual ~BlockForwardsolver() {;}

      /** you guessed it, solve the system  */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const ;
    private:
      TSFArray<TSFLinearSolver> solvers_;
    };

}

#endif
