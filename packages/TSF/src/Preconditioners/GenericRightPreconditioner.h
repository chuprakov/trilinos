#ifndef GENERICRIGHTPRECONDITIONER_H
#define GENERICRIGHTPRECONDITIONER_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFPreconditionerBase.h"

namespace TSF
{
  using std::string;

  /** \ingroup Preconditioner
   * Implements a right preconditioner.
   */

  class GenericRightPreconditioner : public TSFPreconditionerBase
    {
    public:
      /** empty ctor */
      GenericRightPreconditioner(const TSFLinearOperator& op) : op_(op) {;}
      /** TUVD */
      virtual ~GenericRightPreconditioner(){;}

      /** Right preconditioner */
      virtual TSFLinearOperator right() const {return op_;}

      /** return true if this preconditioner has a nontrivial right component */
      virtual bool hasRight() const {return true;}

      /** print to a string */
      virtual string toString() const {return "RightPrecond[" + op_.toString() + "]" ;}

    private:
      TSFLinearOperator op_;
    };


}

#endif
