/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFVECTORSPACE_HPP
#define TSFVECTORSPACE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFCoreVectorSpace.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   *
   */
  template <class Scalar>
  class VectorSpace : public Handle<const TSFCore::VectorSpace<Scalar> >,
                      public TSFCore::VectorSpace<Scalar>
  {
  public:
    /** Empty ctor */
    VectorSpace();

    /** Construct with a raw pointer to a concrete type. Note that the
     * type must implement the Handleable interface. */
    VectorSpace(Handleable<const TSFCore::VectorSpace<Scalar> >* rawPtr);

    /** Construct with a smart pointer to a TSFCore vector space */
    VectorSpace(const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& smartPtr);
    
    /** Create a new element of this vector space */
    RefCountPtr<TSFCore::Vector<Scalar> > createMember() const 
    {return ptr()->createMember();}

    /** Return the dimension of the space */
    int dim() const {return ptr()->dim();}

    /** Check compatibility with another space. Implementation note: 
     * we don't know if the argument vec space is a handle to another
     * vector space or the contents of a handle, and we want the operation
     * to work the same in either case. We can make this work as
     * follows: have the argument check compatibility with the contents
     * of this handle. If the argument is a handle, the process 
     * will be repeated, interchanging places again so that both handles
     * are dereferenced. If the argument is not a handle, then it
     * ends up comparing to the concrete contents of this handle, giving the
     * same results. */
    bool isCompatible(const TSFCore::VectorSpace<Scalar>& vecSpc) const 
    {return vecSpc.isCompatible(*ptr());}
  };
  

  template <class Scalar> inline 
  VectorSpace<Scalar>::VectorSpace()
    : Handle<const TSFCore::VectorSpace<Scalar> >()
  {}

  template <class Scalar> inline 
  VectorSpace<Scalar>::VectorSpace(Handleable<const TSFCore::VectorSpace<Scalar> >* ptr)
    : Handle<const TSFCore::VectorSpace<Scalar> >(ptr)
  {}

  template <class Scalar> inline 
  VectorSpace<Scalar>::VectorSpace(const RefCountPtr<const TSFCore::VectorSpace<Scalar> >& smartPtr)
    : Handle<const TSFCore::VectorSpace<Scalar> >(smartPtr)
  {}
}


#endif
