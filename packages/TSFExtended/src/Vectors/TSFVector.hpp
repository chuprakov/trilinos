/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFVECTOR_HPP
#define TSFVECTOR_HPP

#include "TSFConfigDefs.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace TSFExtended
{
  template <class Scalar> class DeferredLinearCombination;


  /** \ingroup Vector
   * Vector:  User-level handle class for vectors.
   *
   * <b> Constructing a Vector </b>
   *
   * Ordinarily, you will never construct a Vector directly from a
   * derived type.  Rather, the createMember() method of
   * VectorSpace is used to build a vector of the appropriate type,
   * for example,
   * \code 
   * Vector x = mySpace.createMember(); 
   * \endcode
   * This hides from you all the ugly details of creating a particular
   * concrete type.
   *
   * You will frequently create an empty vector to be filled in later, 
   * for example,
   * \code
   * Vector y;
   * \endcode
   * Note that this vector isn't just empty, it's null. Not only does it have no
   * values assigned, it does not have a concrete type. An attempt to set an element
   * of or operate on a null vector will result in an error. What you can do is
   * assign another vector to it,
   * \code
   * Vector y;
   * y = x.deepCopy();
   * \endcode
   * or fill it by passing it as a return-by-reference argument to a mathematical
   * operation,
   * \code
   * Vector y;
   * x.scalarMult(3.14159, y);
   * \endcode
   *
   *
   * Another way a Vector can be created is as the result of an operation between
   * existing vectors, for example
   * \code
   * Vector z = a*x + b*y;
   * \endcode
   * See the section on operator overloading for more information on this.
   *
   *
   * <b> Element access </b>
   *
   * There are a number of methods to get or set the value of an element or group
   * of elements at a particular index or group of indices. In a distributed environment,
   * there is a distinction between an element's <i> local </i> index and its
   * <i> global </i> index. All the element access methods of Vector use global
   * indices.
   *
   * Element access methods must make virtual function calls, so they
   * should be used as sparingly as possible. In particular, you
   * should <i> never </i> do mathematical operations on the whole
   * vector using the element access methods. Use one of the
   * predefined operations (e.g., eMult, update) or a custom reduction
   * and transformation operator instead. If you are needing to set
   * elements, such as when creating the load vector in a
   * finite-element problem, it is more efficient to set groups of
   * elements rather than single elements.
   *
   * <b> Overloaded operators </b>
   *
   * Vector supports the full set of overloaded operators for
   * vector and scalar-vector arithmetic. Overloaded operators are not
   * often used in scientific computing because a naive implementation
   * incurs a large performance penalty due to the creation and
   * copying of temporary vectors. Vector operators use a design
   * which avoids the creation of such temporary vectors by flattening
   * the expression tree. Small temporary objects are created in the
   * process, so for small problems (under \f$N \approx 1000 \f$)
   * there is still a noticeable performance penatly.  However, for
   * larger problems the performance penalty is less than a few
   * percent compared to FORTRAN BLAS, and decreases quickly with
   * problem size. */
  template <class Scalar>
    class Vector
    {
    public:
      /** \name Constructors, Destructors, and Assignment Operators */
      //@{
      /** empty ctor. Constructs a null vector */
      Vector();

      /** Construct with an existing smart pointer */
      Vector(RefCountPtr<TSFCore::Vector<Scalar> >& smartPtr);

      /** construct with the result of an operation */
      Vector(const DeferredLinearCombination<Scalar>& op);

      /** construct with the result of an operation */
      Vector(const DeferredCopy<Scalar>& copy);

      /** assign the result of an operation into this vector */
      Vector& operator=(const DeferredLinearCombination<Scalar>& op);

      /** assign the result of an operation into this vector */
      Vector& operator=(const DeferredCopy<Scalar>& copy);
      //@}

      /** \name Vector space information */
      //@{
      /** return the vector space in which this vector lives. */
      const VectorSpace<Scalar>& space() const ;
      //@}

    private:
      RefCountPtr<TSFCore::Vector<Scalar> > ptr_;
    };

  /** \relates Vector
   * write an Vector to an output stream
   */
  inline template <class Scalar> ostream& operator<<(ostream& os, const Vector<Scalar>& x)
    {
      x.print(os);
      return os;
    }

}


#endif
