/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFVECTOR_HPP
#define TSFVECTOR_HPP

#include "TSFConfigDefs.hpp"
#include "TSFHandle.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFLoadableVector.hpp"
#include "TSFAccessibleVector.hpp"
#include "TSFCoreVectorStdOps.hpp"
#include "Teuchos_TimeMonitor.hpp"

namespace TSFExtendedOps
{
  template <class Scalar, class Node1, class Node2> class LCN;
  template <class Scalar> class LC1;
}

namespace TSFExtended
{
  using TSFCore::Index;

  

  /** 
   * User-level vector class. 
   *
   * <h2> Creating vectors </h2>
   *
   * Ordinarily, you will never construct a Vector directly
   * from a derived type.  Rather, the createMember() method of
   * VectorSpace is used to build a vector of the appropriate
   * type, for example,
   * \code 
   * VectorType<double> vecType = new EpetraVectorType();
   * int dimension = 100;
   * VectorSpace<double> space = vecType.createSpace(dimension);
   * Vector<double> x = space.createMember(); 
   * Vector<double> y = space.createMember(); 
   * \endcode 
   * This hides from you all the ugly
   * details of creating a particular concrete type.
   *
   * You will frequently create an empty vector to be filled in later, 
   * for example,
   * \code
   * Vector<double> y;
   * \endcode
   * Note that this vector isn't just empty, it's null. Not only does 
   * it have no values assigned, it does not have a concrete type. An 
   * call a method on a null vector will result in an error. What you 
   * <it>can</it> do with a null vector is
   * <ul>
   * <li> assign another vector to it
   * \code
   * Vector<double> x = space.createVector();
   * Vector<Scalar> y;
   * y = x.copy();
   * \endcode
   * <li> assign the result of a vector operation to it
   * \code
   * Vector<Scalar> z = a*x + b*y;
   * \endcode
   */
  template <class Scalar>
  class Vector : public Handle<TSFCore::Vector<Scalar> >
    {
    public:
      /** \name Constructors, Destructors, and Assignment Operators */
      //@{
      /** empty ctor. Constructs a null vector */
      Vector();

      /** Construct with an existing smart pointer */
      Vector(const RefCountPtr<TSFCore::Vector<Scalar> >& smartPtr);

#ifndef DOXYGEN_DEVELOPER_ONLY
      /** Assign a one-term linear combination 
       * (i.e., a scalar times a vector) 
       * to this vector */
      Vector& operator=(const TSFExtendedOps::LC1<Scalar>& x);

      /** Assign a linear combination of vectors to this vector */
      template<class Node1, class Node2>
      Vector& operator=(const TSFExtendedOps::LCN<Scalar, Node1, Node2>& x);
#endif
      //@}

      /** */
      RefCountPtr<const TSFCore::VectorSpace<Scalar> > space() const 
      {return ptr()->space();}

      /** \name Math operations */
      //@{
      /** Multiply this vector by a constant scalar factor 
       * \code
       * this = alpha * this;
       * \endcode
      */
      Vector<Scalar>& scale(const Scalar& alpha);

      /** 
       * Add a scaled vector to this vector:
       * \code
       * this = this + alpha*x 
       * \endcode
       */
      Vector<Scalar>& update(const Scalar& alpha, const Vector<Scalar>& x);

      /** 
       * Copy the values of another vector into this vector
       * \code
       * this = x
       * \endcode
       */
      Vector<Scalar>& acceptCopyOf(const Vector<Scalar>& x);

      /** 
       * Create a new vector that is a copy of this vector 
       */
      Vector<Scalar> copy() const ;

      /** 
       * Take dot product with another vector
       */
      Scalar dot() const ;
      //@}

      /** \name Element loading interface */
      //@{
      /** set a single element at the given global index */
      void setElement(Index globalIndex, const Scalar& value) 
      {castToLoadable()->setElement(globalIndex, value);}

      /** add to the existing value of 
       * a single element at the given global index */
      void addToElement(Index globalIndex, const Scalar& value) 
      {castToLoadable()->addToElement(globalIndex, value);}

      /** set a group of elements */
      void setElements(size_t numElems, const Index* globalIndices, 
                               const Scalar* values) 
      {castToLoadable()->setElements(numElems, globalIndices, values);}

      /** add to a group of elements */
      void addToElements(size_t numElems, const Index* globalIndices, 
                         const Scalar* values)
      {castToLoadable()->addToElements(numElems, globalIndices, values);}

      /** Do whatever finalization steps are needed by the implementation,
       for instance, synchronizing border elements. The default implementation
      * is a no-op. */
      void finalizeAssembly() {castToLoadable()->finalizeAssembly();}
      //@}

      /** \name Element access interface */
      //@{
      /** get the element at the given global index */
      const Scalar& getElement(Index globalIndex) const 
      {return castToAccessible()->getElement(globalIndex);}
      //@}
      
      /** Get a stopwtach for timing vector operations */
      static RefCountPtr<Time>& opTimer()
      {
        static RefCountPtr<Time> rtn 
          = TimeMonitor::getNewTimer("Low-level vector operations");
        return rtn;
      }
    private:
#ifndef DOXYGEN_DEVELOPER_ONLY
      /** Cross-cast vector pointer to an accessible vector */
      const AccessibleVector<Scalar>* castToAccessible() const ;

      /** Cross-cast vector to a loadable vector */
      LoadableVector<Scalar>* castToLoadable()  ;

      
      
#endif
    };

  template <class Scalar> inline 
  Vector<Scalar>::Vector()
    : Handle<TSFCore::Vector<Scalar> >()
  {}


  template <class Scalar> inline 
  Vector<Scalar>::Vector(const RefCountPtr<TSFCore::Vector<Scalar> >& smartPtr)
    : Handle<TSFCore::Vector<Scalar> >(smartPtr)
  {}

  template <class Scalar> inline 
  const AccessibleVector<Scalar>* Vector<Scalar>::castToAccessible() const
  {
    const AccessibleVector<Scalar>* av 
      = dynamic_cast<const AccessibleVector<Scalar>*>(ptr().get());
    TEST_FOR_EXCEPTION(av==0, std::runtime_error,
                       "Attempted to cast non-accessible vector "
                       << *this << " to an AccessibleVector");
    return av;
  }

  template <class Scalar> inline 
  LoadableVector<Scalar>* Vector<Scalar>::castToLoadable()
  {
    LoadableVector<Scalar>* lv 
      = dynamic_cast<LoadableVector<Scalar>*>(ptr().get());
    TEST_FOR_EXCEPTION(lv==0, std::runtime_error,
                       "Attempted to cast non-loadable vector "
                       << *this << " to a LoadableVector");
    return lv;
  }

  template <class Scalar> inline 
  Vector<Scalar>& Vector<Scalar>::scale(const Scalar& alpha)
  {
    TSFCore::Vector<Scalar>* p = ptr().get();
    {
      TimeMonitor t(*opTimer());
      TSFCore::Vt_S(p, alpha);
    }
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, const Vector<Scalar>& x)
  {
    TSFCore::Vector<Scalar>* p = ptr().get();
    const TSFCore::Vector<Scalar>* px = x.ptr().get();
    {
      TimeMonitor t(*opTimer());
      TSFCore::Vp_StV(p, alpha, *px);
    }
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
  {
    TSFCore::Vector<Scalar>* p = ptr().get();
    const TSFCore::Vector<Scalar>* px = x.ptr().get();
    {
      TimeMonitor t(*opTimer());
      TSFCore::assign(p, *px);
    }
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar> Vector<Scalar>::copy() const 
  {
    Vector<Scalar> rtn = space()->createMember();
    {
      TimeMonitor t(*opTimer());
      rtn.acceptCopyOf(*this);
    }
    return rtn;
  }

  

  
}


#endif
