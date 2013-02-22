// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP
#define STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP

#include "KokkosArray_Cuda.hpp"
#include "Cuda/KokkosArray_Cuda_Parallel.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_BlockCrsMatrix.hpp"
#include "Stokhos_StochasticProductTensor.hpp"

namespace Stokhos {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< typename TensorScalar ,
          typename MatrixScalar ,
          typename VectorScalar ,
          typename TensorType >
class Multiply<
  BlockCrsMatrix<
    StochasticProductTensor<TensorScalar, TensorType, KokkosArray::Cuda>,
    MatrixScalar, KokkosArray::Cuda> ,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>,
  KokkosArray::View<VectorScalar**, KokkosArray::LayoutLeft, KokkosArray::Cuda>
  >
{
public:

  typedef KokkosArray::Cuda                    device_type ;
  typedef device_type::size_type  size_type ;

  typedef StochasticProductTensor< TensorScalar , TensorType , KokkosArray::Cuda > tensor_type ;
  typedef BlockCrsMatrix< tensor_type, MatrixScalar, device_type > matrix_type ;
  typedef KokkosArray::View< VectorScalar** , KokkosArray::LayoutLeft , KokkosArray::Cuda >           vector_type ;


  static void apply( const matrix_type & A ,
                     const vector_type & x ,
                     const vector_type & y )
  {
    typedef BlockCrsMatrix< TensorType, MatrixScalar , KokkosArray::Cuda > base_matrix_type ;

    typedef Multiply< base_matrix_type , vector_type , vector_type >
      base_multiply_type ;

    base_matrix_type base_matrix ;

    base_matrix.values = A.values ;
    base_matrix.graph  = A.graph ;
    base_matrix.block  = A.block.tensor();

    base_multiply_type::apply( base_matrix , x , y );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_CUDA_STOCHASTICPRODUCTTENSOR_HPP */

