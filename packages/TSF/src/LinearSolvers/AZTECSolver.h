#ifndef AZTECSOLVER_H
#define AZTECSOLVER_H

#include "TSFDefs.h"
#include "TSFLinearSolverBase.h"
#include "TSFTimeMonitor.h"
#include "TSFHashtable.h"



#include "AztecOO.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup ConcreteLinearSolvers
   * Aztec solver
   */

  class AZTECSolver : public TSFLinearSolverBase
    {
    public:
      /** construct a AZTEC solver with default options and parameters. */
      AZTECSolver();

      /** construct a AZTEC solver with the given options and parameters. */
      AZTECSolver(const TSFHashtable<int, int>& aztecOptions,
                  const TSFHashtable<int, double>& aztecParameters);

      /** TUVD */
      virtual ~AZTECSolver();

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
      /** Aztec options */
      mutable TSFArray<int> options_;

      /** Aztec parameters */
      mutable TSFArray<double> parameters_;

    };

}


#endif
