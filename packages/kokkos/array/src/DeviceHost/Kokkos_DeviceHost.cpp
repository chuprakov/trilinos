/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#include <Kokkos_DeviceHost.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace {

class HostMemoryImpl {
public:
  Impl::MemoryInfoSet m_allocations ;

  ~HostMemoryImpl();

  static HostMemoryImpl & singleton();
};

HostMemoryImpl & HostMemoryImpl::singleton()
{
  static HostMemoryImpl self ;
  return self ;
}

HostMemoryImpl::~HostMemoryImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "Kokkos::DeviceHost memory leaks:" << std::endl ;
    m_allocations.print( std::cerr );
  }
}

}

/*--------------------------------------------------------------------------*/

void * HostMemory::allocate_memory(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  Impl::MemoryInfo tmp ;

  tmp.m_type  = & type ;
  tmp.m_label = label ;
  tmp.m_size  = member_size ;
  tmp.m_count = member_count ;
  tmp.m_ptr   = calloc( member_size , member_count );

  const bool ok_alloc  = 0 != tmp.m_ptr ;
  const bool ok_insert = ok_alloc && s.m_allocations.insert( tmp );

  if ( ! ok_alloc || ! ok_insert ) {
    std::ostringstream msg ;
    msg << "Kokkos::HostMemory::allocate_memory( " << label
        << " , " << type.name()
        << " , " << member_size
        << " , " << member_count
        << " ) FAILED " ;
    if ( ok_alloc ) { msg << "memory allocation" ; }
    else            { msg << "with internal error" ; }
    throw std::runtime_error( msg.str() );
  }

  return tmp.m_ptr ;
}

void HostMemory::deallocate_memory( void * ptr )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  if ( ! s.m_allocations.erase( ptr ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::HostMemory::deallocate_memory( " << ptr
        << " ) FAILED memory allocated by this device" ;
    throw std::runtime_error( msg.str() );
  }

  free( ptr );
}

void HostMemory::print_memory_view( std::ostream & o )
{
  HostMemoryImpl & s = HostMemoryImpl::singleton();

  s.m_allocations.print( o );
}

/*--------------------------------------------------------------------------*/

unsigned int HostMemory::m_launching_kernel = false ;

void HostMemory::set_dispatch_functor()
{
  if ( m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::HostMemory::set_dispatch_functor FAILED: " );
    msg.append( "kernel dispatch is already in progress, " );
    msg.append( "a recursive call or forgotten 'clear_dispatch_kernel" );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = true ;
}

void HostMemory::clear_dispatch_functor()
{
  if ( ! m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::HostMemory::clear_dispatch_functor FAILED: " );
    msg.append( "no kernel dispatch in progress." );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = false ;
}

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

