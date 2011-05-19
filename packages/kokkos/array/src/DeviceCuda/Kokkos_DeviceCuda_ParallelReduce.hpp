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

#ifndef KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP
#define KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP

#include <DeviceCuda/Kokkos_DeviceCuda_ParallelDriver.hpp>

#include <impl/Kokkos_DeviceCuda_macros.hpp>

#if defined( KOKKOS_MACRO_DEVICE_FUNCTION )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Must have consistent '__shared__' statement across all device kernels.
// Since there may be more than one kernel in a file then have to make this
// a simple array of words.

namespace Kokkos {

template< class FunctorType , class FinalizeType >
class ParallelReduce< FunctorType , FinalizeType , DeviceCuda > {
public:
  typedef ParallelReduce< FunctorType , FinalizeType , DeviceCuda > self_type ;

  typedef DeviceCuda                        device_type ;
  typedef FunctorType                       functor_type ;
  typedef FinalizeType                      finalize_type ;
  typedef device_type::size_type            size_type ;
  typedef typename functor_type::value_type value_type ;

  enum { WarpSize   = Impl::DeviceCudaTraits::WarpSize };
  enum { WarpStride = WarpSize + 1 };

  enum { SharedMemoryBanks = Impl::DeviceCudaTraits::SharedMemoryBanks_13 };

  enum { ValueWordCount = ( sizeof(value_type) + sizeof(size_type) - 1 )
                          / sizeof(size_type) };

  /** \brief  If the reduction value occupies an
   *          exact multiple of shared memory banks
   *          then it must be padded to avoid bank conflicts.
   */
  enum { ValueWordStride = ValueWordCount +
          ( ValueWordCount % SharedMemoryBanks ? 0 : 2 ) };

  //----------------------------------------------------------------------

  const FunctorType  m_work_functor ;
  const FinalizeType m_work_finalize ;
  const size_type    m_work_count ;
  size_type        * m_scratch_space ;
  size_type        * m_scratch_flag ;
  
  //----------------------------------------------------------------------

  static inline
  __device__
  size_type shared_data_offset()
  { return ValueWordStride * ( threadIdx.x + WarpStride * threadIdx.y ); }

  static inline
  __device__
  size_type shared_data_offset( size_type x , size_type y )
  { return ValueWordStride * ( x + WarpStride * y ); }

  static inline
  __device__
  size_type shared_flag_offset()
  { return ValueWordStride * ( WarpStride * blockDim.y - 1 ); }

  //----------------------------------------------------------------------

  static 
  __device__
  void reduce_shared( const size_type used_warp_count )
  {
    typedef volatile value_type * vvp ;

    enum { HalfWarpSize = WarpSize >> 1 };

    extern __shared__ size_type shared_data[];

    // threadIdx.x == index within warp [ 0 .. WarpSize - 1 ]
    // threadIdx.y == which warp        [ 0 .. used_warp_count - 1 ]
  
    // Phase A: Reduce within my warp:
    //          Warp's reads occur before joins and writes
    //          so there is no race condition.
    //          Declare shared data to be volatile to
    //          prevent compiler from introducing a race condition.
    //
    if ( threadIdx.y < used_warp_count && threadIdx.x < HalfWarpSize ) {
      enum { n1  = ValueWordStride * 1 };
      enum { n2  = ValueWordStride * 2 };
      enum { n4  = ValueWordStride * 4 };
      enum { n8  = ValueWordStride * 8 };
      enum { n16 = ValueWordStride * 16 };

      size_type * const data = shared_data + shared_data_offset();
  
      functor_type::join( *((vvp) data), *((vvp)( data + n16 )) );
      functor_type::join( *((vvp) data), *((vvp)( data +  n8 )) );
      functor_type::join( *((vvp) data), *((vvp)( data +  n4 )) );
      functor_type::join( *((vvp) data), *((vvp)( data +  n2 )) );
      functor_type::join( *((vvp) data), *((vvp)( data +  n1 )) );
    }

    // Phase B: Use a single warp to reduce results from each warp.
    //          This requires: used_warp_count <= WarpSize
    //

    __syncthreads();

    if ( 0 == threadIdx.y && threadIdx.x + 1 < used_warp_count ) {
      enum { n1  = WarpStride * ValueWordStride * 1 };
      enum { n2  = WarpStride * ValueWordStride * 2 };
      enum { n4  = WarpStride * ValueWordStride * 4 };
      enum { n8  = WarpStride * ValueWordStride * 8 };
      enum { n16 = WarpStride * ValueWordStride * 16 };

      size_type * const data = shared_data + shared_data_offset( 0 , threadIdx.x );

      if ( threadIdx.x + 2 < used_warp_count ) {
        if ( threadIdx.x + 4 < used_warp_count ) {
          if ( threadIdx.x + 8 < used_warp_count ) {
            if ( threadIdx.x + 16 < used_warp_count ) {
              functor_type::join( *((vvp) data) , *((vvp)( data + n16 )) );
            }
            functor_type::join( *((vvp) data) , *((vvp)( data + n8 )) );
          }
          functor_type::join( *((vvp) data) , *((vvp)( data + n4 )) );
        }
        functor_type::join( *((vvp) data) , *((vvp)( data + n2 )) );
      }
      functor_type::join( *((vvp) data) , *((vvp)( data + n1 )) );
    }
  }

  //----------------------------------------------------------------------

  inline
  __device__
  void reduce_global() const
  {
    enum { WarpIndexMask  = Impl::DeviceCudaTraits::WarpIndexMask };
    enum { WarpIndexShift = Impl::DeviceCudaTraits::WarpIndexShift };

    extern __shared__ size_type shared_data[];

    const size_type thread_id = threadIdx.x + blockDim.x * threadIdx.y ;

    // Phase A: Output block's results to global memory
    //          and then input results into last block.
    //
    // ** REQUIRED: gridDim.x <= blockDim.x * blockDim.y **

    // Coalesced global memory write, wait for write to complete
    // Write by blockIdx.x and
    // Read  by blockIdx.x == ( threadIdx.x + blockDim.x * threadIdx.y )
    // Determine correct location into scratch memory.

    if ( thread_id < ValueWordCount ) {
      const size_type thread_stride = blockDim.x * blockDim.y ;

      size_type * const scratch = m_scratch_space +
        shared_data_offset( blockIdx.x &  WarpIndexMask  /* for threadIdx.x */ ,
                            blockIdx.x >> WarpIndexShift /* for threadIdx.y */ );

      for ( size_type i = thread_id ; i < ValueWordCount ; i += thread_stride ) {
        scratch[i] = shared_data[i] ;
      }

      __threadfence(); // Wait for write to complete
    }

    // Check if this is the last block to finish:

    if ( 0 == thread_id ) {
      // atomicInc returns value prior to increment.

      shared_data[ shared_flag_offset() ] =
        gridDim.x == 1 + atomicInc( m_scratch_flag , gridDim.x + 1 );
    }
    __syncthreads();

    // The last block to finish reduces the results from all blocks.

    if ( shared_data[ shared_flag_offset() ] ) {

      // Let the last thread of the last warp reinitialize the flag
      // for the next reduce operation.
      // This warp is least likely to be doing any more work.

      if ( blockDim.y == threadIdx.y + 1 &&
           blockDim.x == threadIdx.x + 1 ) {
        *(m_scratch_flag) = 0 ;
      }

      // Number of warps required for reduction of per-block contributions.
      const size_type scratch_warp = ( gridDim.x >> WarpIndexShift ) +
                                     ( gridDim.x &  WarpIndexMask ? 1 : 0 );

      // Each warp does a coalesced read of its own data.

      if ( threadIdx.y < scratch_warp ) {
        const size_type upper_bound = gridDim.x + scratch_warp - 1 ;

        // Coalesced global memory read for this warp's data.

        size_type i = shared_data_offset( 0 , threadIdx.y );
        size_type j = i + WarpSize ;

        if ( upper_bound < j ) {
          j = upper_bound ;

          // Only partial data will be read by this warp
          // so initialize the values before reading.
          functor_type::init( *((value_type *)( shared_data + shared_data_offset() )) );
        }

        i = i * ValueWordStride + threadIdx.x ;
        j = j * ValueWordStride ;

        for ( ; i < j ; i += WarpSize ) {
          shared_data[i] = m_scratch_space[i] ;
        }
      }

      // Phase B: Reduce these contributions
      reduce_shared( scratch_warp );
    }
  }

  //----------------------------------------------------------------------

  __device__
  void run_on_device() const
  {
    extern __shared__ device_type::size_type shared_data[];

    value_type & value =
      *( (value_type *)( shared_data + shared_data_offset() ) );

    functor_type::init( value );

    // Phase 1: Reduce to per-thread contributions
    {
      const size_type work_stride = blockDim.x * blockDim.y * gridDim.x ;

      size_type iwork =
        threadIdx.x + blockDim.x * ( threadIdx.y + blockDim.y * blockIdx.x );

      for ( ; iwork < m_work_count ; iwork += work_stride ) {
        m_work_functor( iwork , value );
      }
    }

    // Phase 2: Reduce this block's thread's contributions
    //          to a single reduction value.
    reduce_shared( blockDim.y );

    // Phase 3: Reduce contributions from multiple blocks

    int last_block = 1 == gridDim.x ;

    if ( ! last_block ) {

      reduce_global();

      last_block = shared_data[ shared_flag_offset() ];
    }

    // Phase 4: Thread #0 of last block performs serial finalization
    if ( last_block && 0 == threadIdx.x && 0 == threadIdx.y ) {
      m_work_finalize( value );
    }
  }

  //----------------------------------------------------------------------

  ParallelReduce( const size_type work_count ,
                  const FunctorType & functor ,
                  const FinalizeType & finalize )
    : m_work_functor( functor )
    , m_work_finalize( finalize )
    , m_work_count( work_count )
    , m_scratch_space( device_type::reduce_multiblock_scratch_space() )
    , m_scratch_flag(  device_type::reduce_multiblock_scratch_flag() )
    {}

  //----------------------------------------------------------------------

  static void run( const size_type      work_count ,
                   const FunctorType  & work_functor ,
                   const FinalizeType & work_finalize )
  {
    // Make a copy just like other devices will have to.

    device_type::set_dispatch_functor();

    const self_type tmp( work_count , work_functor , work_finalize );

    device_type::clear_dispatch_functor();

    dim3 block( Impl::DeviceCudaTraits::WarpSize , 
                device_type::reduce_maximum_warp_count( ValueWordStride ) , 1 );

    dim3 grid( 1 , 1 , 1 );

    if ( work_count <= WarpSize * block.y ) { // Need at most one block
      while ( work_count <= WarpSize * ( block.y >> 1 ) ) { block.y >>= 1 ; }
    }
    else {
      // At most one block per final reduction thread count.
      grid.x = WarpSize * block.y ;

      // Reduce to actual number of blocks needed
      while ( work_count <= WarpSize * block.y * ( grid.x >> 1 ) )
      { grid.x >>= 1 ; }
    }

    const size_type shmem_size =
      sizeof(size_type) * ( ValueWordStride * ( WarpStride * block.y - 1 ) + 1 );

    Impl::device_cuda_run<self_type>( tmp , block , grid , shmem_size );
  }
};

} // namespace Kokkos

#endif /* defined( KOKKOS_MACRO_DEVICE_FUNCTION ) */
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif /* KOKKOS_DEVICECUDA_PARALLELREDUCE_HPP */

