/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <Kokkos_Cuda.hpp>
#include <Kokkos_CudaSpace.hpp>

#include <Cuda/Kokkos_Cuda_Internal.hpp>
#include <impl/Kokkos_MemoryTracking.hpp>
#include <impl/Kokkos_Error.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

DeepCopy<CudaSpace,CudaSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

DeepCopy<HostSpace,CudaSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

DeepCopy<CudaSpace,HostSpace>
  ::DeepCopy( void * dst , const void * src , size_t n )
{ CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) ); }

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
namespace {

class CudaMemoryTracking {
public:

  enum SpaceTag { CudaSpaceTag , CudaUVMSpaceTag , CudaHostPinnedSpaceTag };

  struct Attribute {

    Kokkos::Impl::cuda_texture_object_type m_tex_obj ;

    Attribute() : m_tex_obj(0) {}

    ~Attribute()
      {
        if ( m_tex_obj ) {
          cudaDestroyTextureObject( m_tex_obj );
          m_tex_obj = 0 ;
        }
      }

    cudaError create( void * const                  arg_alloc_ptr
                    , size_t const                  arg_byte_size
                    , cudaChannelFormatDesc const & arg_desc
                    )
    {
      cudaError cuda_status = cudaSuccess ;

      if ( 0 == m_tex_obj ) {
 
        cuda_status = cudaDeviceSynchronize();

        struct cudaResourceDesc resDesc ;
        struct cudaTextureDesc  texDesc ;

        memset( & resDesc , 0 , sizeof(resDesc) );
        memset( & texDesc , 0 , sizeof(texDesc) );

        resDesc.resType                = cudaResourceTypeLinear ;
        resDesc.res.linear.desc        = arg_desc ;
        resDesc.res.linear.sizeInBytes = arg_byte_size ;
        resDesc.res.linear.devPtr      = arg_alloc_ptr ;

        cuda_status = cudaCreateTextureObject( & m_tex_obj , & resDesc, & texDesc, NULL);

        if ( cudaSuccess == cuda_status ) { cuda_status = cudaDeviceSynchronize(); }

        if ( ( cudaSuccess == cuda_status ) && ( 0 == m_tex_obj ) ) {
          throw std::logic_error(std::string("FAILED assumption that Cuda texture objects are non-zero"));
        }
      }

      return cuda_status ;
    }
  };

  typedef          Kokkos::Impl::MemoryTracking< Attribute >         tracking_type ;
  typedef typename Kokkos::Impl::MemoryTracking< Attribute >::Entry  entry_type ;

private:

  tracking_type   m_tracking ;
  SpaceTag const  m_space_tag ;


  cudaError cuda_malloc( void ** ptr , size_t byte_size ) const
    {
      cudaError result = cudaSuccess ;

      switch( m_space_tag ) {
      case CudaSpaceTag :
        result = cudaMalloc( ptr , byte_size );
        break ;
      case CudaUVMSpaceTag :
#if defined( CUDA_VERSION ) && ( 6000 <= CUDA_VERSION )
        result = cudaMallocManaged( ptr, byte_size, cudaMemAttachGlobal );
#else
        Kokkos::Impl::throw_runtime_exception( std::string("CUDA VERSION does not support UVM") );
#endif
        break ;
      case CudaHostPinnedSpaceTag :
        result = cudaHostAlloc( ptr , byte_size , cudaHostAllocDefault );
        break ;
      }

      return result ;
    }

  cudaError cuda_free( void * ptr ) const
    {
      cudaError result = cudaSuccess ;

      switch( m_space_tag ) {
      case CudaSpaceTag :
      case CudaUVMSpaceTag :
        result = cudaFree( ptr );
        break ;
      case CudaHostPinnedSpaceTag :
        result = cudaFreeHost( ptr );
        break ;
      }
      return result ;
    }

public :

  CudaMemoryTracking( const SpaceTag arg_tag , const char * const arg_label )
    : m_tracking(  arg_label )
    , m_space_tag( arg_tag )
    {}

  void print( std::ostream & oss , const std::string & lead ) const
    { m_tracking.print( oss , lead ); }

  const char * query_label( const void * ptr ) const
    {
      static const char error[] = "<NOT FOUND>" ;
      entry_type * const entry = m_tracking.query( ptr );
      return entry ? entry->label() : error ;
    }

  void * allocate(
    const std::string    & label ,
    const std::type_info & scalar_type ,
    const size_t           scalar_size ,
    const size_t           scalar_count )
  {
    void * ptr = 0 ;

    const size_t byte_size = scalar_size * scalar_count ;

    if ( byte_size ) {

      const bool ok_parallel = ! HostSpace::in_parallel();

      cudaError cuda_status = cudaSuccess ;

      if ( ok_parallel ) {

        cuda_status = cudaDeviceSynchronize();

        if ( cudaSuccess == cuda_status ) { cuda_status = CudaMemoryTracking::cuda_malloc( & ptr , byte_size ); }
        if ( cudaSuccess == cuda_status ) { cuda_status = cudaDeviceSynchronize(); }
      }

      if ( ok_parallel && ( cudaSuccess == cuda_status ) ) {
        m_tracking.insert( label , ptr , scalar_size , scalar_count );
      }
      else {
        std::ostringstream msg ;
        msg << m_tracking.label()
            << "::allocate( "
            << label
            << " , " << scalar_type.name()
            << " , " << scalar_size
            << " , " << scalar_count
            << " ) FAILURE : " ;
        if ( ! ok_parallel ) {
          msg << "called within a parallel functor" ;
        }
        else {
          msg << " CUDA ERROR \"" << cudaGetErrorString(cuda_status) << "\"" ;
        }
        Kokkos::Impl::throw_runtime_exception( msg.str() );
      }
    }

    return ptr ;
  }

  void decrement( const void * ptr )
  {
    const bool ok_parallel = ! HostSpace::in_parallel();

    cudaError cuda_status = cudaSuccess ;

    if ( ok_parallel ) {

      cuda_status = cudaDeviceSynchronize();

      void * const alloc_ptr = ( cudaSuccess == cuda_status ) ? m_tracking.decrement( ptr ) : (void *) 0 ;

      if ( alloc_ptr ) {
        if ( cudaSuccess == cuda_status ) { cuda_status = CudaMemoryTracking::cuda_free( alloc_ptr ); }
        if ( cudaSuccess == cuda_status ) { cuda_status = cudaDeviceSynchronize(); }
      }
    }

    if ( ( ! ok_parallel ) || ( cudaSuccess != cuda_status ) ) {
      std::ostringstream msg ;
      msg << m_tracking.label() << "::decrement( " << ptr << " ) FAILURE : " ;
      if ( ! ok_parallel ) {
        msg << "called within a parallel functor" ;
      }
      else {
        msg << " CUDA ERROR \"" << cudaGetErrorString(cuda_status) << "\"" ;
      }
      std::cerr << msg.str() << std::endl ;
    }
  }

  void increment( const void * ptr )
    {
      const bool ok_parallel = ! HostSpace::in_parallel();

      if ( ok_parallel ) {
        m_tracking.increment( ptr );
      }
      else {
        std::ostringstream msg ;
        msg << m_tracking.label() << "::increment(" << ptr
            << ") FAILURE :called within a parallel functor" ;
        Kokkos::Impl::throw_runtime_exception( msg.str() );
      }
    }


  inline
  void texture_object_attach( const void * const            arg_ptr
                            , const cudaChannelFormatDesc & arg_desc
                            , ::cudaTextureObject_t * const arg_tex_obj
                            , void const           ** const arg_alloc_ptr
                            , int                   * const arg_offset
                            )
    {
      static const size_t max_array_len = 1 << 28 ;

      *arg_tex_obj   = 0 ;
      *arg_alloc_ptr = 0 ;
      *arg_offset    = 0 ;

      if ( arg_ptr ) {

        // Can only create texture object on device architure 3.0 or better
        const bool ok_dev_arch = 300 <= Cuda::device_arch();
        const bool ok_parallel = ok_dev_arch && ! HostSpace::in_parallel();

        entry_type * const entry = ok_parallel ? m_tracking.query( arg_ptr ) : (entry_type *) 0 ;

        const size_t offset = entry ? ( reinterpret_cast<const char*>(arg_ptr) -
                                        reinterpret_cast<const char*>(entry->m_alloc_ptr) ) : 0 ;

        const bool ok_offset = entry     && ( 0 == ( offset % entry->m_type_size ) );
        const bool ok_count  = ok_offset && ( entry->m_array_len < max_array_len );

        cudaError cuda_status = cudaSuccess ;

        if ( ok_count ) {
          cuda_status = entry->m_attribute.create( entry->m_alloc_ptr
                                                 , entry->m_type_size * entry->m_array_len
                                                 , arg_desc );
        }

        if ( ok_count && ( cudaSuccess == cuda_status ) ) {
          *arg_tex_obj   = entry->m_attribute.m_tex_obj ;
          *arg_alloc_ptr = entry->m_alloc_ptr ;
          *arg_offset    = offset / entry->m_type_size ;
        }
        else {
          std::ostringstream msg ;
          msg << m_tracking.label()
              << "::cuda_texture_object_attach(" << arg_ptr << ") FAILED :" ;
          if ( ! ok_dev_arch ) {
            msg << " cuda architecture " << Cuda::device_arch()
                << " does not support texture objects" ;
          }
          else if ( ! ok_parallel ) {
            msg << " called within a parallel functor" ;
          }
          else if ( 0 == entry ) {
            msg << " pointer not tracked" ;
          }
          else if ( ! ok_offset ) {
            msg << " pointer not properly aligned" ;
          }
          else if ( ! ok_count ) {
            msg << " array too large for texture object" ;
          }
          else {
            msg << " CUDA ERROR \"" << cudaGetErrorString(cuda_status) << "\"" ;
          }
          Kokkos::Impl::throw_runtime_exception( msg.str() );
        }
      }
    }
};

//----------------------------------------------------------------------------

CudaMemoryTracking &
cuda_space_singleton()
{
#if defined( KOKKOS_USE_CUDA_UVM )
  static CudaMemoryTracking s( CudaMemoryTracking::CudaUVMSpaceTag , "Kokkos::CudaSpace");
#else
  static CudaMemoryTracking s( CudaMemoryTracking::CudaSpaceTag , "Kokkos::CudaSpace");
#endif
  return s ;
}

CudaMemoryTracking &
cuda_uvm_space_singleton()
{
  static CudaMemoryTracking s( CudaMemoryTracking::CudaUVMSpaceTag , "Kokkos::CudaUVMSpace");
  return s ;
}

CudaMemoryTracking &
cuda_host_pinned_space_singleton()
{
  static CudaMemoryTracking s( CudaMemoryTracking::CudaHostPinnedSpaceTag , "Kokkos::CudaHostPinnedSpace");
  return s ;
}

}
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

void * CudaSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  return Impl::cuda_space_singleton().allocate( label , scalar_type , scalar_size , scalar_count );
}

void CudaSpace::decrement( const void * ptr )
{
  Impl::cuda_space_singleton().decrement( ptr );
}


void CudaSpace::increment( const void * ptr )
{
  Impl::cuda_space_singleton().increment( ptr );
}

void CudaSpace::print_memory_view( std::ostream & oss )
{
  Impl::cuda_space_singleton().print( oss , std::string("  ") );
}

std::string CudaSpace::query_label( const void * p )
{
  return std::string( Impl::cuda_space_singleton().query_label(p) );
}

void CudaSpace::texture_object_attach( const void * const            arg_ptr
                                     , ::cudaChannelFormatDesc const & arg_desc
                                     , ::cudaTextureObject_t * const arg_tex_obj
                                     , void const           ** const arg_alloc_ptr
                                     , int                   * const arg_offset
                                     )
{
  Impl::cuda_space_singleton().texture_object_attach( arg_ptr , arg_desc , arg_tex_obj , arg_alloc_ptr , arg_offset );
}

void CudaSpace::access_error()
{
  const std::string msg("Kokkos::CudaSpace::access_error attempt to execute Cuda function from non-Cuda space" );

  Kokkos::Impl::throw_runtime_exception( msg );
}

void CudaSpace::access_error( const void * const ptr )
{
  std::ostringstream msg ;
  msg << "Kokkos::CudaSpace::access_error:" ;
  msg << " attempt to access Cuda-data labeled(" ;
  msg << query_label( ptr ) ;
  msg << ") from non-Cuda execution" ;
  Kokkos::Impl::throw_runtime_exception( msg.str() );
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

void * CudaUVMSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  return Impl::cuda_uvm_space_singleton().allocate( label , scalar_type , scalar_size , scalar_count );
}

void CudaUVMSpace::decrement( const void * ptr )
{
  Impl::cuda_uvm_space_singleton().decrement( ptr );
}


void CudaUVMSpace::increment( const void * ptr )
{
  Impl::cuda_uvm_space_singleton().increment( ptr );
}

void CudaUVMSpace::print_memory_view( std::ostream & oss )
{
  Impl::cuda_uvm_space_singleton().print( oss , std::string("  ") );
}

std::string CudaUVMSpace::query_label( const void * p )
{
  return std::string( Impl::cuda_uvm_space_singleton().query_label(p) );
}

void CudaUVMSpace::texture_object_attach( const void * const            arg_ptr
                                        , ::cudaChannelFormatDesc const & arg_desc
                                        , ::cudaTextureObject_t * const arg_tex_obj
                                        , void const           ** const arg_alloc_ptr
                                        , int                   * const arg_offset
                                        )
{
  Impl::cuda_uvm_space_singleton().texture_object_attach( arg_ptr , arg_desc , arg_tex_obj , arg_alloc_ptr , arg_offset );
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

void * CudaHostPinnedSpace::allocate(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count )
{
  return Impl::cuda_host_pinned_space_singleton().allocate( label , scalar_type , scalar_size , scalar_count );
}

void CudaHostPinnedSpace::decrement( const void * ptr )
{
  Impl::cuda_host_pinned_space_singleton().decrement( ptr );
}


void CudaHostPinnedSpace::increment( const void * ptr )
{
  Impl::cuda_host_pinned_space_singleton().increment( ptr );
}

void CudaHostPinnedSpace::print_memory_view( std::ostream & oss )
{
  Impl::cuda_host_pinned_space_singleton().print( oss , std::string("  ") );
}

std::string CudaHostPinnedSpace::query_label( const void * p )
{
  return std::string( Impl::cuda_host_pinned_space_singleton().query_label(p) );
}

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

