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

#ifndef TSFVECTORTYPE_HPP
#define TSFVECTORTYPE_HPP

#include "TSFHandle.hpp"
#include "TSFVectorTypeExtensions.hpp"
#include "TSFVectorSpaceDecl.hpp"

namespace TSFExtended
{
  using namespace Teuchos;
  /**
   * Vector type objects are used by the application code to create
   * vector spaces and operators of a given type.
   */
  template <class Scalar>
  class VectorType : public Handle<VectorTypeExtensions<Scalar> >
  {
  public:
    HANDLE_CTORS(VectorType<Scalar>, VectorTypeExtensions<Scalar>);
   

    /** Create a vector space in which all elements are replicated on
     * all processors. This is used when creating multivector-based
     * operators. */
    VectorSpace<Scalar> createReplicatedSpace(int dimension) const ;

    /** create a distributed vector space.
     * @param dimension the dimension of the space
     * @param nLocal number of indices owned by the local processor
     * @param locallyOwnedIndices array of indices owned by this processor  
     */
    VectorSpace<Scalar> createSpace(int dimension, 
                                    int nLocal,
                                    const int* locallyOwnedIndices) const ;

    /**
     * Create an empty matrix of type compatible with this vector type,
     * sized according to the given domain and range spaces.
     */
    LinearOperator<Scalar>
    createMatrix(const VectorSpace<Scalar>& domain,
                 const VectorSpace<Scalar>& range) const ;
                                                      
    
  };




  template <class Scalar> inline 
  VectorSpace<Scalar> VectorType<Scalar>::createReplicatedSpace(int dimension) const
  {
    const TSFCore::VectorSpaceFactory<Scalar>* f
      = dynamic_cast<const TSFCore::VectorSpaceFactory<Scalar>*>(ptr().get());
    
    TEST_FOR_EXCEPTION(f==0, runtime_error, 
                       "failed cast to TSFCore::VectorSpaceFactory in "
                       "VectorType<Scalar>::createReplicatedSpace()");
    return f->createVecSpc(dimension);
  }

  template <class Scalar> inline 
  VectorSpace<Scalar> VectorType<Scalar>::createSpace(int dimension,
                                                      int nLocal,
                                                      const int* locallyOwnedIndices) const
  {
    return ptr()->createSpace(dimension, nLocal, locallyOwnedIndices);
  }

  template <class Scalar> inline
  LinearOperator<Scalar>
  VectorType<Scalar>::createMatrix(const VectorSpace<Scalar>& domain,
                                   const VectorSpace<Scalar>& range) const
  {
    return ptr()->createMatrix(domain, range);
  }
}

#endif
