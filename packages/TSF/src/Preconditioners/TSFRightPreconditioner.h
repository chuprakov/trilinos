#ifndef TSFRIGHTPRECONDITIONER_H
#define TSFRIGHTPRECONDITIONER_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFPreconditionerBase.h"

namespace TSF
{
  using std::string;

  // VEH
  /** \ingroup Preconditioner
   * Stores a linear operator M1^-1 that will be applied as
   * a right preconditioner.
   */

  class TSFRightPreconditioner : public TSFPreconditionerBase
    {
    public:
      /** empty ctor */
      TSFRightPreconditioner(const TSFLinearOperator& rightOp);
      /** TUVD */
      virtual ~TSFRightPreconditioner(){;}

      /** Right preconditioner */
      virtual TSFLinearOperator right() const {return rightOp_;}

      /** tell the world that we have a right preconditioner */
      virtual bool hasRight() const {return true;}

      /** print to a string */
      virtual string toString() const ;
    private:
      TSFLinearOperator rightOp_;
    };


}

#endif
