#ifndef TSFVECTORTYPEEXTENSIONS_HPP
#define TSFVECTORTYPEEXTENSIONS_HPP

#include "TSFHandle.hpp"
#include "TSFCoreVectorSpaceFactory.hpp"
#include "TSFVectorSpace.hpp"
#include "TSFLinearOperator.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * VectorTypeExtensions provides extensions to 
   * the TSFCore::VectorSpaceFactory
   * object, appropriate to the interface with applications codes.
   * TSFCore::VectorSpaceFactory has a method to create a small serial
   * vector space for use in multivector operators, however, that is 
   * insufficient for the needs of applications.
   *
   * <h4> Notes for subclass developers </h4>
   *
   * Subclasses should also derive from some subclass of
   * TSFCore::VectorSpaceFactory in order to use the createReplicatedSpace()
   * method. 
   *
   * Because applications may use the clean syntax
   * \code
   * VectorType vt = new MyVectorType();
   * \endcode
   * subclasses should also derive from
   * Handleable<VectorTypeBase<Scalar> and implement the getRcp() method.
   *
   * Subclasses might optionally implement the Printable and Describable
   * interfaces
   *
   */
  template <class Scalar>
  class VectorTypeExtensions 
  {
  public:

    /** create a distributed vector space.
     * @param dimension the dimension of the space 
     * @param nLocal number of indices owned by the local processor
     * @param locallyOwnedIndices array of indices owned by this processor  
     */
    virtual RefCountPtr<const TSFCore::VectorSpace<Scalar> >
    createSpace(int dimension, 
                int nLocal,
                const int* locallyOwnedIndices) const = 0 ;

    
    /**
     * Create an empty matrix of type compatible with this vector type,
     * sized according to the given domain and range spaces.
     */
    virtual LinearOperator<Scalar>
    createMatrix(const VectorSpace<Scalar>& domain,
                 const VectorSpace<Scalar>& range) const = 0 ;
  };
  
}

#endif
