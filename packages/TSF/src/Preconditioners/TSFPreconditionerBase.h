#ifndef TSFPRECONDITIONERBASE_H
#define TSFPRECONDITIONERBASE_H

#include "TSFDefs.h"
#include "TSFVector.h"
#include "TSFLinearOperator.h"

namespace TSF
{
  using std::string;

  /** \ingroup Preconditioner
   * Base class for preconditioners. A general preconditioner object
   * is split into a left preconditioner M1^-1 and a right
   * preconditioner M2^-1. To solve A x = b, we define the auxiliary
   * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y.
   * Having y, we can quickly recover x by applying M2^-1 to y.
   *
   * The base class implements neither a left nor a right preconditioner,
   * and throws an error upon a call to left() or right(). Derived classes
   * should override one of more of these.
   */

  class TSFPreconditionerBase
    {
    public:
      /** empty ctor */
      TSFPreconditionerBase(){;}
      /** TUVC */
      virtual ~TSFPreconditionerBase(){;}

      /** Left preconditioner */
      virtual TSFLinearOperator left() const ;

      /** Right preconditioner */
      virtual TSFLinearOperator right() const ;

      /** return true if this preconditioner has a nontrivial left component */
      virtual bool hasLeft() const {return false;}

      /** return true if this preconditioner has
       * a nontrivial right component */
      virtual bool hasRight() const {return false;}

      /** print to a string */
      virtual string toString() const = 0 ;
    private:
    };


}

#endif
