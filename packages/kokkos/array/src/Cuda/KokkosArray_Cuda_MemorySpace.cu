/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <KokkosArray_Cuda.hpp>
#include <Cuda/KokkosArray_Cuda_MemorySpace.hpp>
#include <Cuda/KokkosArray_Cuda_Internal.hpp>
#include <impl/KokkosArray_MemoryTracking.hpp>

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

namespace {

class CudaMemoryImpl {
public:
  MemoryTracking m_allocations ;

  CudaMemoryImpl();
  ~CudaMemoryImpl();

  static CudaMemoryImpl & singleton();
};

CudaMemoryImpl::CudaMemoryImpl()
  : m_allocations()
{}

CudaMemoryImpl & CudaMemoryImpl::singleton()
{
  static CudaMemoryImpl self ;
  return self ;
}

CudaMemoryImpl::~CudaMemoryImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "KokkosArray::Cuda memory leaks:" << std::endl ;
    m_allocations.print( std::cerr , std::string("  ") );
  }
}

}

/*--------------------------------------------------------------------------*/

void CudaMemorySpace::copy_to_device_from_device(
  void * dst, void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

void CudaMemorySpace::copy_to_device_from_host(
  void * dst, void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}

void CudaMemorySpace::copy_to_host_from_device(
  void * dst, void * src, size_t n )
{
  CUDA_SAFE_CALL( cudaMemcpy( dst , src , n , cudaMemcpyDefault ) );
}



/*--------------------------------------------------------------------------*/

void * CudaMemorySpace::allocate(
  const std::string    & label ,
  const std::type_info & value_type ,
  const size_t           value_size ,
  const size_t           value_count )
{
  CudaMemoryImpl & s = CudaMemoryImpl::singleton();

  const size_t size = value_size * value_count ;

  void * ptr = 0 ;

  if ( value_size * value_count ) {
    bool ok = true ;

    if ( ok ) ok = cudaSuccess == cudaMalloc( & ptr , size );
    if ( ok ) ok = 0 != ptr ;
    if ( ok ) ok = cudaSuccess == cudaMemset( ptr , 0 , size );
    if ( ok ) ok = cudaSuccess == cudaThreadSynchronize();

    if ( ! ok ) {
      std::ostringstream msg ;
      msg << "KokkosArray::Impl::CudaMemorySpace::allocate( "
          << label
          << " , " << value_type.name()
          << " , " << value_size
          << " , " << value_count
          << " ) FAILED memory allocation" ;
      throw std::runtime_error( msg.str() );
    }

    s.m_allocations.track( ptr, & value_type, value_size, value_count, label );
  }

  return ptr ;
}

#if ! defined( __CUDA_ARCH__ )

void CudaMemorySpace::increment( const void * ptr )
{
  if ( 0 != ptr && m_memory_view_tracking ) {
    CudaMemoryImpl & s = CudaMemoryImpl::singleton();

    s.m_allocations.increment( ptr );
  }
}

void CudaMemorySpace::decrement( const void * ptr )
{
  if ( 0 != ptr && m_memory_view_tracking ) {

    CudaMemoryImpl & s = CudaMemoryImpl::singleton();

    MemoryTracking::Info info = s.m_allocations.decrement( ptr );

    if ( 0 == info.count ) {
      const bool failed =
        cudaSuccess != cudaFree( const_cast<void*>( info.begin ) );

      if ( failed ) {
        std::string msg("KokkosArray::Impl::CudaMemorySpace::decrement() failed cudaFree");
        throw std::runtime_error( msg );
      }
    }
  }
}

#endif

void CudaMemorySpace::print_memory_view( std::ostream & o )
{
  CudaMemoryImpl & s = CudaMemoryImpl::singleton();

  s.m_allocations.print( o , std::string("  ") );
}


size_t CudaMemorySpace::preferred_alignment(
  size_t value_size , size_t value_count )
{
  const size_t align = ( Impl::CudaTraits::WarpSize * value_size ) / sizeof(Cuda::size_type);
  const size_t rem = value_count % align ;
  if ( rem ) value_count += align - rem ;
  return value_count ;
}

/*--------------------------------------------------------------------------*/

int CudaMemorySpace::m_memory_view_tracking = true ;

void CudaMemorySpace::disable_memory_view_tracking()
{
  if ( ! m_memory_view_tracking ) {
    std::string msg ;
    msg.append( "KokkosArray::Impl::CudaMemory::disable_memory_view_tracking ");
    msg.append( " FAILED memory_view_tracking already disabled" );
    throw std::runtime_error( msg );
  }
  m_memory_view_tracking = false ;
}

void CudaMemorySpace::enable_memory_view_tracking()
{
  if ( m_memory_view_tracking ) {
    std::string msg ;
    msg.append( "KokkosArray::Impl::CudaMemory::enable_memory_view_tracking ");
    msg.append( " FAILED memory_view_tracking already enabled" );
    throw std::runtime_error( msg );
  }
  m_memory_view_tracking = true ;
}

/*--------------------------------------------------------------------------*/

} // namespace Impl
} // namespace KokkosArray

