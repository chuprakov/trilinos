/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#ifndef KOKKOS_DEVICECUDA_PARALLELFOR_HPP
#define KOKKOS_DEVICECUDA_PARALLELFOR_HPP

namespace Kokkos {
namespace {

template< class ParallelForDriver >
__global__
void cuda_parallel_for_driver()
{
  ParallelForDriver::run_on_device();
}


};

template< class FunctorType >
class ParallelFor< FunctorType , DeviceCuda > {
public:
  typedef DeviceCuda             device_type ;
  typedef device_type::size_type size_type ;
  typedef ParallelFor< FunctorTYpe , DeviceCuda > self_type ;

  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ParallelFor( const size_type work_count , const FunctorType & functor )
    : m_work_functor( functor )
    , m_work_count( work_count )
    {}

  __device__
  void run_on_device()
  {
    const self_type * const self = functor();
    const size_type work_stride = blockDim.x * blockDim.y * gridDim.x ;

    size_type iwork = threadIdx.x + blockDim.x * (
                      threadIdx.y + blockDim.y * blockIdx.x );

    for ( ; iwork < self->m_work_count ; iwork += work_stride ) {
      self->m_functor( iwork );
    }
  }


  static void run( const size_type     work_count ,
                   const FunctorType & work_functor )
  {
    // Make a copy just like other devices will have to.

    device_type::set_dispatch_functor();

    const self_type tmp( work_count , work_functor );

    device_type::clear_dispatch_functor();

    dim3 block( DeviceCuda::Traits::WarpSize , 
                device_type::parallel_for_warp_count() , 1 );

    dim3 grid( device_type::maximum_grid_count() , 1 , 1 );

    // Reduce grid count until just enough blocks for the work.

    while ( work_count <= block.x * block.y * ( grid.x >> 1 ) ) {
      grid.x >>= 1 ; }
    }

    CudaDevice::load_constant_buffer( this , sizeof(self_type) );

    cuda_parallel_for_driver< self_type > <<< grid , block >>>();
  }
};

} // namespace Kokkos

#endif /* KOKKOS_DEVICECUDA_PARALLELFOR_HPP */

