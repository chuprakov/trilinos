#ifndef TSFVECTORTYPE_HPP
#define TSFVECTORTYPE_HPP

#include "TSFHandle.hpp"
#include "TSFVectorTypeBase.hpp"
#include "TSFVectorSpace.hpp"

namespace TSFExtended
{
  using namespace Teuchos;
  /**
   * User-level handle class for vector types. 
   */
  template <class Scalar>
  class VectorType : public Handle<TSFCore::VectorSpaceFactory<Scalar> >
  {
  public:

    /** */
    VectorSpace<Scalar> createSpace(int dimension) const ;

    /** create a vector space in which the local processor owns
     * indices \f$[firstLocal, firstLocal+nLocal]f$. */
    VectorSpace<Scalar> createSpace(int dimension, 
                                    int nLocal,
                                    int firstLocal) const ;
    
    /** create a vector space in which the given local indices are owned by 
     * this processor */
    VectorSpace<Scalar> createSpace(int dimension, 
                                     int nLocal,
                                     const int* localIndices) const ;
    
    /** create a vector space in which the given local indices are owned
     * by this processor, and the given ghost indices are available but
     * not owned. */
    VectorSpace<Scalar> createSpace(int dimension, 
                                      int nLocal,
                                      const int* localIndices,
                                      int nGhost,
                                      const int* ghostIndices) const ;
  };
  
}

#endif
