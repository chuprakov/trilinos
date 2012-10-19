/*
//@HEADER
// ************************************************************************
//
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSARRAY_CUDA_PARALLEL_HPP
#define KOKKOSARRAY_CUDA_PARALLEL_HPP

//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

struct CudaTraits {
  enum { WarpSize       = 32      /* 0x0020 */ };
  enum { WarpIndexMask  = 0x001f  /* Mask for warpindex */ };
  enum { WarpIndexShift = 5       /* WarpSize == 1 << WarpShift */ };
  enum { SharedMemoryBanks_13 = 16 /* Compute device 1.3 */ };
  enum { SharedMemoryBanks_20 = 32 /* Compute device 2.0 */ };
  enum { UpperBoundGridCount = 65535 /* Hard upper bound */ };
  enum { ConstantMemoryCapacity = 0x010000 /* 64k bytes */ };
  enum { ConstantMemoryCache    = 0x002000 /*  8k bytes */ };

  typedef unsigned long
    ConstantGlobalBufferType[ ConstantMemoryCapacity / sizeof(unsigned long) ];

  enum { ConstantMemoryUseThreshold = 0x000100 /* 256 bytes */ };

  static inline
#if defined( __CUDACC__ )
  __device__ __host__
#endif
  CudaSpace::size_type warp_count( CudaSpace::size_type i )
    { return ( i + WarpIndexMask ) >> WarpIndexShift ; }

  static inline
#if defined( __CUDACC__ )
  __device__ __host__
#endif
  CudaSpace::size_type warp_align( CudaSpace::size_type i )
    {
      enum { Mask = ~CudaSpace::size_type( WarpIndexMask ) };
      return ( i + WarpIndexMask ) & Mask ;
    }
};

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

namespace KokkosArray {
namespace Impl {

CudaSpace::size_type cuda_internal_maximum_warp_count();
CudaSpace::size_type cuda_internal_maximum_grid_count();
CudaSpace::size_type cuda_internal_maximum_shared_words();

CudaSpace::size_type * cuda_internal_scratch_flags( CudaSpace::size_type size );
CudaSpace::size_type * cuda_internal_scratch_space( CudaSpace::size_type size );
CudaSpace::size_type * cuda_internal_scratch_unified( CudaSpace::size_type size );

template< typename ValueType >
inline
__device__
void cuda_internal_atomic_add( ValueType & update , ValueType input )
{ atomicAdd( & update , input ); }

inline
__device__
void cuda_internal_atomic_add( double & update , double input )
{
  typedef unsigned long long int UInt64 ;

  UInt64 * const address = reinterpret_cast<UInt64*>( & update );
  UInt64 test ;
  union UType { double d ; UInt64 i ; } value ;

  value.i = *address ; // Read existing value

  do {
    test = value.i ;
    value.d += input ;
    value.i = atomicCAS( address , test , value.i );
  } while ( value.i != test );
}

} // namespace Impl
} // namespace KokkosArray

/** \brief  Access to constant memory on the device */
__device__ __constant__
KokkosArray::Impl::CudaTraits::ConstantGlobalBufferType
kokkos_impl_cuda_constant_memory_buffer ;

template< typename T >
inline
__device__
T * kokkos_impl_cuda_shared_memory()
{ extern __shared__ KokkosArray::CudaSpace::size_type sh[]; return (T*) sh ; }

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
// See section B.17 of Cuda C Programming Guide Version 3.2
// for discussion of
//   __launch_bounds__(maxThreadsPerBlock,minBlocksPerMultiprocessor)
// function qualifier which could be used to improve performance.
//----------------------------------------------------------------------------
// Maximize L1 cache and minimize shared memory:
//   cudaFuncSetCacheConfig(MyKernel, cudaFuncCachePreferL1 );
// For 2.0 capability: 48 KB L1 and 16 KB shared
//----------------------------------------------------------------------------

template< class DriverType >
__global__
static void cuda_parallel_launch_constant_memory()
{
  const DriverType & driver =
    *((const DriverType *) kokkos_impl_cuda_constant_memory_buffer );

  driver();
}

template< class DriverType >
__global__
static void cuda_parallel_launch_local_memory( const DriverType driver )
{
  driver();
}

template < class DriverType ,
           bool Large = ( CudaTraits::ConstantMemoryUseThreshold < sizeof(DriverType) ) >
struct CudaParallelLaunch ;

template < class DriverType >
struct CudaParallelLaunch< DriverType , true > {

  inline
  CudaParallelLaunch( const DriverType & driver ,
                      const dim3       & grid ,
                      const dim3       & block ,
                      const int          shmem )
  {
    // Copy functor to constant memory on the device
    cudaMemcpyToSymbol( kokkos_impl_cuda_constant_memory_buffer , & driver , sizeof(DriverType) );

    // Invoke the driver function on the device
    cuda_parallel_launch_constant_memory< DriverType ><<< grid , block , shmem >>>();
  }
};

template < class DriverType >
struct CudaParallelLaunch< DriverType , false > {

  inline
  CudaParallelLaunch( const DriverType & driver ,
                      const dim3       & grid ,
                      const dim3       & block ,
                      const int          shmem )
  {
    cuda_parallel_launch_local_memory< DriverType ><<< grid , block , shmem >>>( driver );
  }
};

//----------------------------------------------------------------------------

template< typename DstType , typename SrcType  >
class CudaParallelCopy ;

template< typename Type >
class CudaParallelCopy<Type,Type> {
public:
  CudaParallelCopy( Type * dst , const Type * src , CudaSpace::size_type count )
  {
    CUDA_SAFE_CALL( cudaMemcpy( dst , src , count * sizeof(Type) ,
                                cudaMemcpyDefault ) );
  }
};

template< typename DstType , typename SrcType  >
class CudaParallelCopy {
public:

        DstType * const m_dst ;
  const SrcType * const m_src ;
  const CudaSpace::size_type m_count ;
        CudaSpace::size_type m_stride ;

  inline
  __device__
  void operator()(void) const
  {
    CudaSpace::size_type i = threadIdx.x + blockDim.x * blockIdx.x ;
    for ( ; i < m_count ; i += m_stride ) {
      m_dst[i] = (DstType) m_src[i] ;
    }
  }

  CudaParallelCopy( DstType * dst , const SrcType * src ,
                    CudaSpace::size_type count )
    : m_dst( dst ), m_src( src ), m_count( count )
    {
      const CudaSpace::size_type grid_max = cuda_internal_maximum_grid_count();

      const dim3 block( CudaTraits::WarpSize * cuda_internal_maximum_warp_count(), 1, 1);

      dim3 grid( ( ( count + block.x - 1 ) / block.x ) , 1 , 1 );

      if ( grid_max < grid.x ) grid.x = grid_max ;

      m_stride = grid.x * block.x ;

      cuda_parallel_launch_local_memory< CudaParallelCopy ><<< grid , block >>>( *this );
    }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------

#endif /* defined( __CUDACC__ ) */

#endif /* #define KOKKOSARRAY_CUDA_PARALLEL_HPP */

