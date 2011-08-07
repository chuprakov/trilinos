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

#ifndef KOKKOS_DEVICETBB_MULTIVECTORVIEW_HPP
#define KOKKOS_DEVICETBB_MULTIVECTORVIEW_HPP

#include <Kokkos_MultiVectorView.hpp>

#include <Kokkos_DeviceTBB_macros.hpp>
#include <impl/Kokkos_MultiVectorView_macros.hpp>
#include <Kokkos_DeviceClear_macros.hpp>

namespace Kokkos {
namespace Impl {

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceTBB , true , DeviceHost , true >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceTBB >  & dst ,
                   const MultiVectorView< ValueType , DeviceHost > & src )
  {
    typedef MultiVectorDeepCopy< ValueType , DeviceTBB , true ,
                                             DeviceTBB , true > functor_type ;

    parallel_for( dst.length() * dst.count() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

template< typename ValueType >
class MultiVectorDeepCopy< ValueType , DeviceHost , true , DeviceTBB , true >
{
public:
  static void run( const MultiVectorView< ValueType , DeviceHost > & dst ,
                   const MultiVectorView< ValueType , DeviceTBB >  & src )
  {
    typedef MultiVectorDeepCopy< ValueType , DeviceTBB , true ,
                                             DeviceTBB , true > functor_type ;

    parallel_for( dst.length() * dst.count() ,
                  functor_type( dst.m_memory.ptr_on_device() ,
                                src.m_memory.ptr_on_device() ) );
  }
};

} // namespace Impl
} // namespace Kokkos

#endif /* #ifndef KOKKOS_DEVICETBB_MULTIVECTORVIEW_HPP */


