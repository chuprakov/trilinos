#ifndef TSFTWOSIDEDPRECONDITIONER_H
#define TSFTWOSIDEDPRECONDITIONER_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFPreconditionerBase.h"

namespace TSF
{
  using std::string;

  /** \ingroup Preconditioner
   * Stores two linear operators, used for the left and right preconditioners.
   */

  class TSFTwoSidedPreconditioner : public TSFPreconditionerBase
    {
    public:
      /** empty ctor */
      TSFTwoSidedPreconditioner(const TSFLinearOperator& leftOp,
                                const TSFLinearOperator& rightOp);
      /** TUVD */
      virtual ~TSFTwoSidedPreconditioner(){;}

      /** Left preconditioner */
      virtual TSFLinearOperator left() const {return leftOp_;}

      /** Right preconditioner */
      virtual TSFLinearOperator right() const {return rightOp_;}

      /** inform the world that we have a nontrivial left precond */
      virtual bool hasLeft() const {return true;}

      /** inform the world that we have a nontrivial right precond */
      virtual bool hasRight() const {return true;}

      /** print to a string */
      virtual string toString() const ;
    private:
      TSFLinearOperator leftOp_;
      TSFLinearOperator rightOp_;
    };
}

#endif
