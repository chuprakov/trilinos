/* @HEADER@ */
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

namespace TSFExtended
{
  using TSFCore::Index;

  template <class Scalar, class Node1, class Node2> class LCN;
  template <class Scalar> class LC1;

  /** 
   * User-level vector class. 
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

      /** */
      Vector& operator=(const LC1<Scalar>& x);

      /** */
      template<class Node1, class Node2>
      Vector& operator=(const LCN<Scalar, Node1, Node2>& x);
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
      
      
    private:
      /** Cross-cast vector pointer to an accessible vector */
      const AccessibleVector<Scalar>* castToAccessible() const ;

      /** Cross-cast vector to a loadable vector */
      LoadableVector<Scalar>* castToLoadable()  ;
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
    TSFCore::Vt_S(p, alpha);
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar>& Vector<Scalar>::update(const Scalar& alpha, const Vector<Scalar>& x)
  {
    TSFCore::Vector<Scalar>* p = ptr().get();
    const TSFCore::Vector<Scalar>* px = x.ptr().get();
    TSFCore::Vp_StV(p, alpha, *px);
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar>& Vector<Scalar>::acceptCopyOf(const Vector<Scalar>& x)
  {
    TSFCore::Vector<Scalar>* p = ptr().get();
    const TSFCore::Vector<Scalar>* px = x.ptr().get();
    TSFCore::assign(p, *px);
    return *this;
  }

  template <class Scalar> inline 
  Vector<Scalar> Vector<Scalar>::copy() const 
  {
    Vector<Scalar> rtn = space()->createMember();
    rtn.acceptCopyOf(*this);
    return rtn;
  }

  

  
}


#endif
