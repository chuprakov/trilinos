#ifndef TSFNONLINEAROPERATOR_H
#define TSFNONLINEAROPERATOR_H

#include "TSFDefs.h"
#include "TSFLinearOperator.h"
#include "TSFDeferredLinearCombination.h"
#include "TSFNonlinearOperatorBase.h"

namespace TSF
{






  /** \ingroup NonlinearOperator
   * User-level handle for nonlinear operator objects.
   *
   */

  class TSFNonlinearOperator
    {
    public:
      /** empty ctor constructs a null nonlinear operator.
       * This is primarily for
       * use with templated container classes. */
      TSFNonlinearOperator();
      /** create a TSFNonlinearOperator from a pointer to a subtype. */
      TSFNonlinearOperator(TSFNonlinearOperatorBase* ptr);
      /** create a TSFNonlinearOperator from a TSFLinearOperator. */
      TSFNonlinearOperator(const TSFLinearOperator& linOp);

      /** return domain space */
      const TSFVectorSpace& domain() const ;
      /** return range space */
      const TSFVectorSpace& range() const ;

      /** apply the operator to a vector. The input vector is tested to ensure
       * that it is in the domain of the operator */
      void apply(const TSFVector& arg, TSFVector& out) const ;

      /** get a linear operator representing the derivative */
      TSFLinearOperator derivative(const TSFVector& evalPt) const ;

      /** scale by a constant */
      TSFNonlinearOperator operator*(const TSFReal& scale) const ;

      /** addition of two operators */
      TSFNonlinearOperator operator+(const TSFNonlinearOperator& op) const ;

      /** subtraction of two operators */
      TSFNonlinearOperator operator-(const TSFNonlinearOperator& op) const ;

      /** negation of an operator */
      TSFNonlinearOperator operator-() const ;

      /** sum of an operator and a constant-valued operator  */
      TSFNonlinearOperator operator+(const TSFVector& v) const ;

      /** difference of an operator and a constant-valued operator  */
      TSFNonlinearOperator operator-(const TSFVector& v) const ;

      /** composition of two operators */
      TSFNonlinearOperator compose(const TSFNonlinearOperator& op) const ;

      /** access to the pointer to the linear operator implementation.
       *  For developer use only. */
      TSFSmartPtr<TSFNonlinearOperatorBase>& getPtr() {return ptr_;}

      /** access to the pointer to the linear operator implementation.
       *  For developer use only. */
      const TSFSmartPtr<TSFNonlinearOperatorBase>& getPtr() const {return ptr_;}

      /** write to a stream */
      void print(ostream& os) const ;
    private:

      /* pointer to a concrete type */
      TSFSmartPtr<TSFNonlinearOperatorBase> ptr_;
    };

  /** \relates TSFNonlinearOperator left-multiplication by a scalar */
  inline TSFNonlinearOperator operator*(const TSFReal& a,
                                        const TSFNonlinearOperator& op)
    {
      return op*a;
    }

  /** \relates TSFNonlinearOperator addition with a vector */
  inline TSFNonlinearOperator operator+(const TSFVector& v,
                                        const TSFNonlinearOperator& op)
    {
      return op+v;
    }

  /** \relates TSFNonlinearOperator subtraction from a vector */
  inline TSFNonlinearOperator operator-(const TSFVector& v,
                                        const TSFNonlinearOperator& op)
    {
      return -op+v;
    }

  /** \relates TSFNonlinearOperator write to a stream */
  inline ostream& operator<<(ostream& os, const  TSFNonlinearOperator& op)
    {
      op.print(os);
      return os;
    }

}

#endif
