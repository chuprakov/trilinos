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

#ifndef KOKKOS_HOSTDEVICEFOR_HPP
#define KOKKOS_HOSTDEVICEFOR_HPP

#include <algorithm>
#include <TPI.h>
#include <Kokkos_ParallelFor.hpp>
#include <Kokkos_HostTPI.hpp>

namespace Kokkos {

template< class FunctorType >
struct ParallelFor< FunctorType , HostTPI > {

  typedef HostTPI::size_type size_type ;

  const FunctorType m_functor ;
  const size_type   m_work_count ;

  static void run_functor_on_tpi( TPI_Work * work )
  {
    const ParallelFor & self = *((const ParallelFor *) work->info );

    const size_type work_inc   = ( self.m_work_count + work->count - 1 ) / work->count ;
    const size_type work_begin = work_inc * work->rank ;
    const size_type work_end   = std::max( work_begin + work_inc , self.m_work_count );

    for ( size_type iwork = work_begin ; iwork < work_end ; ++iwork ) {
      self.m_functor( iwork );
    }
  }

  inline
  ParallelFor( const size_type     arg_work_count ,
               const FunctorType & arg_functor )
    : m_work_count( arg_work_count )
    , m_functor(    arg_functor )
    { TPI_Run_threads( & run_functor_on_tpi , this , 0 ); }

  static void run( const size_type work_count , const FunctorType & functor )
  { ParallelFor( work_count , functor ); }
};

}

#endif /* KOKKOS_HOSTDEVICEFOR_HPP */

