#ifndef TSFPRECONDITIONER_H
#define TSFPRECONDITIONER_H

#include "TSFDefs.h"
#include "TSFPreconditionerBase.h"
#include "TSFSmartPtr.h"
#include "TSFLinearOperator.h"

namespace TSF
{


  /** \ingroup Core
   * User-level interface to preconditioners. A general preconditioner object
   * is split into a left preconditioner M1^-1 and a right
   * preconditioner M2^-1. To solve A x = b, we define the auxiliary
   * system M2^-1 y = x, and solve M1^-1 A M2^-1 y = M1^-1 b to obtain y.
   * Having y, we can quickly recover x by applying M2^-1 to y.
   */

  class TSFPreconditioner
    {
    public:
      /** empty ctor constructs a null linear operator. This is primarily for
       * use with templated container classes. */
      TSFPreconditioner();
      /** create a TSFLinearOperator from a pointer to a subtype. */
      TSFPreconditioner(TSFPreconditionerBase* ptr);

      /** Left preconditioner */
      TSFLinearOperator left() const {return ptr_->left();}

      /** Right preconditioner */
      TSFLinearOperator right() const {return ptr_->right();}

      /** return true if this preconditioner has both left and
       * right components. */
      bool isTwoSided() const {return hasLeft() && hasRight();}

      /** return true if this preconditioner has a nontrivial left component */
      bool hasLeft() const ;

      /** return true if this preconditioner has
       * a nontrivial right component */
      bool hasRight() const ;

      /** return true if this preconditioner has neither left nor
       * right operators defined */
      bool isIdentity() const {return !hasLeft() && !hasRight();}

      /** print as a string */
      string toString() const {return ptr_->toString();}


    private:
      /* pointer to a concrete type */
      TSFSmartPtr<TSFPreconditionerBase> ptr_;
    };

  /** \relates TSFPreconditioner write to stream */
  inline ostream& operator<<(ostream& os, const TSFPreconditioner& p)
    {
      return os << p.toString();
    }

  /** \relates TSFPreconditioner write to string */
  inline string toString(const TSFPreconditioner& p)
    {
      return p.toString();
    }

}

#endif
