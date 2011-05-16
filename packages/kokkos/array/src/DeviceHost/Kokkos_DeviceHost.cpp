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

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

#include <Kokkos_DeviceHost.hpp>
#include <impl/Kokkos_MemoryInfo.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace {

class DeviceHostImpl {
public:
  Impl::MemoryInfoSet m_allocations ;

  ~DeviceHostImpl();

  static DeviceHostImpl & singleton();
};

DeviceHostImpl & DeviceHostImpl::singleton()
{
  static DeviceHostImpl self ;
  return self ;
}

DeviceHostImpl::~DeviceHostImpl()
{
  if ( ! m_allocations.empty() ) {
    std::cerr << "Kokkos::DeviceHost memory leaks:" << std::endl ;
    m_allocations.print( std::cerr );
  }
}

}

/*--------------------------------------------------------------------------*/

void * DeviceHost::allocate_memory(
  const std::string    & label ,
  const std::type_info & type ,
  const size_t member_size ,
  const size_t member_count )
{
  DeviceHostImpl & s = DeviceHostImpl::singleton();

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
    msg << "Kokkos::DeviceHost::allocate_memory( " << label
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

void DeviceHost::deallocate_memory( void * ptr )
{
  DeviceHostImpl & s = DeviceHostImpl::singleton();

  if ( ! s.m_allocations.erase( ptr ) ) {
    std::ostringstream msg ;
    msg << "Kokkos::DeviceHost::deallocate_memory( " << ptr
        << " ) FAILED memory allocated by this device" ;
    throw std::runtime_error( msg.str() );
  }

  free( ptr );
}

void DeviceHost::print_memory_view( std::ostream & o )
{
  DeviceHostImpl & s = DeviceHostImpl::singleton();

  s.m_allocations.print( o );
}

/*--------------------------------------------------------------------------*/

unsigned int DeviceHost::m_launching_kernel = false ;

void DeviceHost::set_dispatch_functor()
{
  if ( m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceHost::set_dispatch_functor FAILED: " );
    msg.append( "kernel dispatch is already in progress, " );
    msg.append( "a recursive call or forgotten 'clear_dispatch_kernel" );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = true ;
}

void DeviceHost::clear_dispatch_functor()
{
  if ( ! m_launching_kernel ) {
    std::string msg ;
    msg.append( "Kokkos::DeviceHost::clear_dispatch_functor FAILED: " );
    msg.append( "no kernel dispatch in progress." );
    throw std::runtime_error( msg );
  }
  m_launching_kernel = false ;
}

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

