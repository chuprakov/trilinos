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

#ifndef KOKKOSARRAY_HOST_PARALLELREDUCE_HPP
#define KOKKOSARRAY_HOST_PARALLELREDUCE_HPP

#include <KokkosArray_ParallelReduce.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>

#include <algorithm>
#include <vector>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< typename ValueType >
class ParallelReduceFunctorValue< ValueType , Host >
{
public:
  typedef ValueType value_type ;

  ParallelReduceFunctorValue() {}

  inline void operator()( const value_type & ) const {}

  value_type result() const
  {
    value_type * const ptr = (value_type*) Host::root_reduce_scratch();
    return *ptr ;
  }
};

template< typename MemberType >
class ParallelReduceFunctorValue< MemberType[] , Host >
{
public:
  typedef MemberType    value_type[] ;
  const HostSpace::size_type value_count ;

  inline void operator()( const MemberType [] ) const {}

  explicit
  ParallelReduceFunctorValue( HostSpace::size_type n )
    : value_count(n)
    {}

  void result( value_type result ) const
  {
    MemberType * const ptr = (MemberType *) Host::root_reduce_scratch();

    for ( HostSpace::size_type i = 0 ; i < value_count ; ++i ) result[i] = ptr[i] ;
  }
};

//----------------------------------------------------------------------------

template< class FunctorType , class ValueOper , class FinalizeType , class WorkSpec >
class ParallelReduce< FunctorType , ValueOper , FinalizeType , Host , WorkSpec >
  : public HostThreadWorker
{
public:

  typedef ReduceOperator< ValueOper , FinalizeType >  reduce_oper ;
  typedef          Host::size_type         size_type ;
  typedef typename ValueOper::value_type  value_type ;

  const FunctorType   m_work_functor ;
  const reduce_oper   m_reduce ;
  const size_type     m_work_count ;

  void execute_on_thread( HostThread & this_thread ) const
  {
#if defined( __INTEL_COMPILER )
    enum { work_align = is_same<WorkSpec,VectorParallel>::value && 
                        power_of_two<HostSpace::WORK_ALIGNMENT>::value
                      ? HostSpace::WORK_ALIGNMENT : 1 };
    enum { work_mask  = work_align - 1 };
#else
    enum { work_align = 1 };
    enum { work_mask  = 0 };
#endif

    // Iterate this thread's work

    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    if ( ! work_mask ) {
      // This thread's reduction value, initialized
      m_reduce.init( this_thread.reduce_data() );

      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_work_functor( iwork , m_reduce.reference( this_thread.reduce_data() ) );
      }
    }
#if defined( __INTEL_COMPILER )
    else {

#pragma simd
#pragma ivdep
      for ( size_type j = 0 ; j < HostSpace::WORK_ALIGNMENT ; ++j ) {
        m_reduce.init( this_thread.reduce_data() , j );
      }

#pragma simd vectorlength(work_align)
#pragma ivdep
      for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
        m_work_functor( iwork , m_reduce.reference( this_thread.reduce_data() , iwork & mask_align ) );
      }

      m_reduce.template join< HostSpace::WORK_ALIGNMENT >( this_thread.reduce_data() );
    }
#endif

    // End the routine with a reduction.
    end_reduce( this_thread , m_reduce );
  }

  ParallelReduce( const size_type      work_count ,
                  const FunctorType  & functor ,
                  const FinalizeType & finalize )
    : m_work_functor( functor )
    , m_reduce( finalize )
    , m_work_count( work_count )
    {
      Host::resize_reduce_scratch( m_reduce.value_size() );
      HostThreadWorker::execute();
      m_reduce.finalize( Host::root_reduce_scratch() );
    }
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class FunctorType , class ReduceOper >
class HostMultiFunctorParallelReduceMember ;

template< class ReduceOper >
struct HostMultiFunctorParallelReduceMember<void,ReduceOper> {

  virtual ~HostMultiFunctorParallelReduceMember() {}

  virtual void apply( HostThread & , const ReduceOper & ) const = 0 ;
};


template< class FunctorType , class ReduceOper >
class HostMultiFunctorParallelReduceMember
  : public HostMultiFunctorParallelReduceMember<void,ReduceOper> {
public:
  typedef Host::size_type size_type ;
    
  const FunctorType m_work_functor ;
  const size_type   m_work_count ;

  ~HostMultiFunctorParallelReduceMember() {}

  HostMultiFunctorParallelReduceMember(
    const FunctorType & work_functor ,
    const size_type work_count )
    : m_work_functor( work_functor )
    , m_work_count(   work_count )
    {}
    
  // virtual method
  void apply( HostThread & this_thread ,
              const ReduceOper & reduce ) const
  {
    const std::pair<size_type,size_type> range =
      this_thread.work_range( m_work_count );

    for ( size_type iwork = range.first ; iwork < range.second ; ++iwork ) {
      m_work_functor( iwork , reduce.reference( this_thread.reduce_data() ) );
    }
  }
};  

} // namespace Impl
  
template< class ValueOper , class FinalizeType >
class MultiFunctorParallelReduce< ValueOper , FinalizeType , Host >
  : public Impl::HostThreadWorker
{
public:

  typedef Impl::ReduceOperator< ValueOper , FinalizeType > reduce_oper ;
  typedef          Host::size_type         size_type ;
  typedef typename ValueOper::value_type  value_type ;
  typedef Impl::HostMultiFunctorParallelReduceMember<void,reduce_oper> worker_type ;

  typedef std::vector< worker_type * > MemberContainer ;

  typedef typename MemberContainer::const_iterator MemberIterator ;

  MemberContainer m_member_functors ;
  reduce_oper     m_reduce ;

  void execute_on_thread( Impl::HostThread & this_thread ) const
  {
    // This thread's reduction value, initialized
    m_reduce.init( this_thread.reduce_data() );

    for ( MemberIterator m  = m_member_functors.begin() ;
                         m != m_member_functors.end() ; ++m ) {
      (*m)->apply( this_thread , m_reduce );
    }

    // End the routine with a reduction
    end_reduce( this_thread , m_reduce );
  }

public:

  MultiFunctorParallelReduce( const FinalizeType & finalize )
    : m_member_functors()
    , m_reduce( finalize )
    { }

  ~MultiFunctorParallelReduce()
  {
    while ( ! m_member_functors.empty() ) {
      delete m_member_functors.back();
      m_member_functors.pop_back();
    }
  }

  template< class FunctorType >
  void push_back( const size_type work_count , const FunctorType & functor )
  {
    typedef Impl::HostMultiFunctorParallelReduceMember<FunctorType,reduce_oper> member_work_type ;

    worker_type * const m = new member_work_type( functor , work_count );

    m_member_functors.push_back( m );
  }

  void execute() const
  {
    Host::resize_reduce_scratch( m_reduce.value_size() );
    Impl::HostThreadWorker::execute();
    m_reduce.finalize( Host::root_reduce_scratch() );
  }
};

} // namespace KokkosArray

#endif /* KOKKOSARRAY_HOST_PARALLELREDUCE_HPP */

