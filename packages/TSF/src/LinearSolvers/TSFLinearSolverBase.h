#ifndef TSFLINEARSOLVERBASE_H
#define TSFLINEARSOLVERBASE_H

#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFParameterList.h"
#include "TSFParameterListImplem.h"

namespace TSF
{
  using std::ostream;

  /** \ingroup LinearSolverSubtypes
   */

  class TSFLinearSolverBase
    {
    public:
      /** */
      TSFLinearSolverBase(const TSFParameterList& params = TSFParameterList("Empty"));
      /** TUVD */
      virtual ~TSFLinearSolverBase();

      /**
       * Solve the system with the given RHS, returning the solution
       * by reference
       * argument. The return value is true if the solve succeeded, false
       * if it failed.
       */
      virtual bool solve(const TSFLinearOperator& op,
                         const TSFVector& rhs,
                         TSFVector& soln) const = 0 ;

      /** */
      void setVerbosityLevel(int v) {verbosity_ = v;}

      /** */
      virtual TSFParameterList defaultParameters() const ;

    protected:
      int verbosity_;

      TSFParameterList params_;
    private:

    };

}


#endif
