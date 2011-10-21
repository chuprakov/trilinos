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

#ifndef KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP
#define KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP

#include <Kokkos_ParallelFor.hpp>

#include <algorithm>
#include <vector>

namespace Kokkos {
namespace Impl {

template< class FunctorType >
class ParallelFor< FunctorType , DevicePthread > : public DevicePthreadWorker {
public:

  typedef DevicePthread::size_type size_type ;

  const FunctorType m_work_functor ;

private:

  virtual
  void execute_on_thread( Impl::DevicePthreadController & this_thread ) const
  {
    // Iterate this thread's work
    size_type iwork = DevicePthreadWorker::m_work_portion * this_thread.rank();

    const size_type work_end =
      std::min( iwork + DevicePthreadWorker::m_work_portion ,
                        DevicePthreadWorker::m_work_count );

    for ( ; iwork < work_end ; ++iwork ) {
      m_work_functor( iwork );
    }

    this_thread.barrier();
  }

  ParallelFor( const size_type work_count ,
               const FunctorType & functor )
    : DevicePthreadWorker( work_count )
    , m_work_functor( functor )
    {}

public:

  static void execute( const size_type     work_count ,
                       const FunctorType & functor )
  {
    DevicePthread::memory_space::set_dispatch_functor();

    ParallelFor driver( work_count , functor );

    DevicePthread::memory_space::clear_dispatch_functor();

    DevicePthread::execute( driver );
  }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_DEVICEPTHREAD_PARALLELFOR_HPP */

