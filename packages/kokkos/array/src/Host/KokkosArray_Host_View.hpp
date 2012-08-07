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

#ifndef KOKKOS_HOST_VIEW_HPP
#define KOKKOS_HOST_VIEW_HPP

#include <KokkosArray_View.hpp>
#include <Host/KokkosArray_Host_Parallel.hpp>

#include <KokkosArray_Host_macros.hpp>
#include <impl/KokkosArray_ViewOperLeft_macros.hpp>
#include <impl/KokkosArray_ViewOperRight_macros.hpp>
#include <impl/KokkosArray_View_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {

template< class DataType , class LayoutType >
void View< DataType , LayoutType , Host >::create(
  const std::string & label ,
  const View< DataType , LayoutType , Host >::shape_type shape )
{
  const size_t count = Impl::allocation_count( shape );

  oper_type::m_shape = shape ;
  oper_type::m_ptr_on_device = (value_type *)
    memory_space::allocate( label ,
                            typeid(value_type) ,
                            sizeof(value_type) ,
                            count );

  Impl::HostParallelFill<value_type>( oper_type::m_ptr_on_device , 0 , count );
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template< class OutputView , class InputView  , unsigned Rank >
struct HostViewRemap ;

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 8 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( Host::size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
    for ( Host::size_type i7 = 0 ; i7 < output.dimension_7() ; ++i7 ) {
      output(i0,i1,i2,i3,i4,i5,i6,i7) = input(i0,i1,i2,i3,i4,i5,i6,i7);
    }}}}}}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 7 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
    for ( Host::size_type i6 = 0 ; i6 < output.dimension_6() ; ++i6 ) {
      output(i0,i1,i2,i3,i4,i5,i6) = input(i0,i1,i2,i3,i4,i5,i6);
    }}}}}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 6 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
    for ( Host::size_type i5 = 0 ; i5 < output.dimension_5() ; ++i5 ) {
      output(i0,i1,i2,i3,i4,i5) = input(i0,i1,i2,i3,i4,i5);
    }}}}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 5 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
    for ( Host::size_type i4 = 0 ; i4 < output.dimension_4() ; ++i4 ) {
      output(i0,i1,i2,i3,i4) = input(i0,i1,i2,i3,i4);
    }}}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 4 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
    for ( Host::size_type i3 = 0 ; i3 < output.dimension_3() ; ++i3 ) {
      output(i0,i1,i2,i3) = input(i0,i1,i2,i3);
    }}}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 3 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
    for ( Host::size_type i2 = 0 ; i2 < output.dimension_2() ; ++i2 ) {
      output(i0,i1,i2) = input(i0,i1,i2);
    }}}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 2 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
    for ( Host::size_type i1 = 0 ; i1 < output.dimension_1() ; ++i1 ) {
      output(i0,i1) = input(i0,i1);
    }}

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 1 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { HostThreadWorker<void>::execute( *this ); }

  void execute_on_thread( HostThread & this_thread ) const
  {
    std::pair<Host::size_type,Host::size_type> range =
      this_thread.work_range( output.dimension_0() );

    for ( Host::size_type i0 = range.first ; i0 < range.second ; ++i0 ) {
      output(i0) = input(i0);
    }

    this_thread.barrier();
  }
};

template< class OutputView , class InputView >
struct HostViewRemap< OutputView , InputView , 0 >
  : public HostThreadWorker<void>
{
  const OutputView output ;
  const InputView  input ;

  HostViewRemap( const OutputView & arg_out , const InputView & arg_in )
    : output( arg_out ), input( arg_in )
    { *arg_out = *arg_in ; }
};

//----------------------------------------------------------------------------
// Deep copy views with either different value types
// or different layouts.

template< class DataTypeDst , class LayoutDst ,
          class DataTypeSrc , class LayoutSrc >
struct ViewDeepCopy< View< DataTypeDst , LayoutDst , Host > ,
                     View< DataTypeSrc , LayoutSrc , Host > ,
                     true_type  /* Same value_type  */ ,
                     false_type /* Different layout_type */ ,
                     true_type  /* Same rank */ >
{
  typedef View< DataTypeDst , LayoutDst , Host > dst_type ;
  typedef View< DataTypeSrc , LayoutSrc , Host > src_type ;

  static inline
  void apply( const dst_type & dst , const src_type & src )
  {
    assert_shapes_equal_dimension( dst.shape() , src.shape() );

    HostViewRemap< dst_type , src_type , dst_type::Rank >( dst , src );
  }
};

template< class DataTypeDst , class LayoutDst ,
          class DataTypeSrc , class LayoutSrc ,
          class SameLayout >
struct ViewDeepCopy< View< DataTypeDst , LayoutDst , Host > ,
                     View< DataTypeSrc , LayoutSrc , Host > ,
                     false_type /* Different value_type  */ ,
                     SameLayout /* Any layout */ ,
                     true_type  /* Same rank */ >
{
  typedef View< DataTypeDst , LayoutDst , Host > dst_type ;
  typedef View< DataTypeSrc , LayoutSrc , Host > src_type ;

  static inline
  void apply( const dst_type & dst , const src_type & src )
  {
    assert_shapes_equal_dimension( dst.shape() , src.shape() );

    HostViewRemap< dst_type , src_type , dst_type::Rank >( dst , src );
  }
};

} // namespace Impl

//----------------------------------------------------------------------------

template< typename ValueType , class LayoutSrc >
inline
void deep_copy( ValueType & dst ,
                const View< ValueType , LayoutSrc , Host > & src )
{
  typedef View< ValueType , LayoutSrc , Host > src_type ;
  typedef typename src_type::shape_type        src_shape ;

  typedef typename Impl::assert_shape_is_rank_zero< src_shape >::type ok_rank ;

  dst = *src ;
}

template< typename ValueType , class LayoutDst >
inline
void deep_copy( const View< ValueType , LayoutDst , Host > & dst ,
                const ValueType & src )
{
  typedef View< ValueType , LayoutDst , Host > dst_type ;
  typedef typename dst_type::shape_type        dst_shape ;

  typedef typename Impl::assert_shape_is_rank_zero< dst_shape >::type ok_rank ;

  *dst = src ;
}

} // namespace KokkosArray

#endif /* #ifndef KOKKOS_HOST_VIEW_HPP */


