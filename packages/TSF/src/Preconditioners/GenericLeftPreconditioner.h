#ifndef GENERICLEFTPRECONDITIONER_H
#define GENERICLEFTPRECONDITIONER_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFLinearOperator.h"

namespace TSF
{
  using std::string;

  /** \ingroup Preconditioner
   * Implements a left preconditioner.
   */

  class GenericLeftPreconditioner : public TSFPreconditionerBase
    {
    public:
      /** empty ctor */
      GenericLeftPreconditioner(const TSFLinearOperator& op) : op_(op) {;}
      /** TUVD */
      virtual ~GenericLeftPreconditioner(){;}

      /** Left preconditioner */
      virtual TSFLinearOperator left() const {return op_;}

      /** return true if this preconditioner has a nontrivial left component */
      virtual bool hasLeft() const {return true;}

      /** print to a string */
      virtual string toString() const {return "LeftPrecond[" + op_.toString() + "]" ;}

    private:
      TSFLinearOperator op_;
    };


}

#endif
