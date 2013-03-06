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

#ifndef KOKKOSARRAY_VIEWRIGHT_HPP
#define KOKKOSARRAY_VIEWRIGHT_HPP

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstViewType >
struct ViewAssignment<
  DstViewType ,
  typename DstViewType::memory_space ,
  typename enable_if< (
    ( is_same< typename DstViewType::array_layout , LayoutRight >::value )
    &&
    ( is_same< typename DstViewType::memory_traits , MemoryManaged >::value )
  ) >::type >
{
  typedef typename DstViewType::shape_type shape_type ;

private:

  typedef typename DstViewType::memory_space  memory_space ;
  typedef typename DstViewType::layout_type   layout_type ;

  static inline
  size_t stride( const shape_type & shape )
  {
    const size_t block_count =
                 shape.N1 * shape.N2 * shape.N3 *
      shape.N4 * shape.N5 * shape.N6 * shape.N7 ;

    return
      memory_space::preferred_alignment( shape.scalar_size , block_count );
  }

  static inline
  void allocate( DstViewType & dst , const std::string & label )
  {
    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    const size_t allocation_count = dst.m_shape.N0 * dst.m_stride ;

    dst.m_ptr_on_device = (typename DstViewType::scalar_type *)
      memory_space::allocate( label ,
                              typeid(typename DstViewType::scalar_type) ,
                              sizeof(typename DstViewType::scalar_type) ,
                              allocation_count );
  }

public:

  // Same data type, same layout, different device; used to create a mirror.
  template< class SrcViewType >
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  typename enable_if<
                    is_same< DstViewType , typename SrcViewType::HostMirror >::value
                  >::type * = 0 )
  {
    dst.m_shape = src.m_shape ;
    dst.m_stride = src.m_stride ;
    allocate( dst , "mirror" );
  }

  ViewAssignment( DstViewType & dst , const std::string & label , const shape_type shape )
  {
    dst.m_shape = shape ;
    dst.m_stride = stride( shape );

    allocate( dst , label );
  }

  ViewAssignment( DstViewType & dst , const std::string & label ,
                  const size_t n0 = 0 ,
                  const size_t n1 = 0 ,
                  const size_t n2 = 0 ,
                  const size_t n3 = 0 ,
                  const size_t n4 = 0 ,
                  const size_t n5 = 0 ,
                  const size_t n6 = 0 ,
                  const size_t n7 = 0 )
  {
    shape_type::assign( dst.m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );

    dst.m_stride = stride( dst.m_shape );

    allocate( dst , label );
  }
};

} /* namespace Impl */
} /* namespace KokkosArray */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 1 )
  ) , unsigned_<1> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 )
  {
    assert_shape_bounds( src.shape() , i0 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i0 ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 2 )
  ) , unsigned_<2> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device = src.m_ptr_on_device + i1 + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 3 )
  ) , unsigned_<3> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i2 + src.m_shape.N2 * (
      i1 ) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 4 )
  ) , unsigned_<4> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i3 + src.m_shape.N3 * (
      i2 + src.m_shape.N2 * (
      i1 )) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 5 )
  ) , unsigned_<5> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i4 + src.m_shape.N4 * (
      i3 + src.m_shape.N3 * (
      i2 + src.m_shape.N2 * (
      i1 ))) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 6 )
  ) , unsigned_<6> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i5 + src.m_shape.N5 * (
      i4 + src.m_shape.N4 * (
      i3 + src.m_shape.N3 * (
      i2 + src.m_shape.N2 * (
      i1 )))) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 7 )
  ) , unsigned_<7> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 , const unsigned i6 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 , i6 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i6 + src.m_shape.N6 * (
      i5 + src.m_shape.N5 * (
      i4 + src.m_shape.N4 * (
      i3 + src.m_shape.N3 * (
      i2 + src.m_shape.N2 * (
      i1 ))))) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

template< class DstViewType , class SrcViewType >
struct ViewAssignment< DstViewType , SrcViewType ,
  typename enable_if< (
    ( ViewAssignment_Compatible<
        typename DstViewType::view_traits ,
        typename SrcViewType::view_traits >::compatible_value )
    &&
    ( DstViewType::Rank == 0 )
    &&
    ( is_same< typename SrcViewType::array_layout , LayoutRight >::value )
    &&
    ( SrcViewType::Rank == 8 )
  ) , unsigned_<8> >::type >
{
  KOKKOSARRAY_INLINE_FUNCTION
  ViewAssignment( DstViewType & dst , const SrcViewType & src ,
                  const unsigned i0 , const unsigned i1 , const unsigned i2 , const unsigned i3 ,
                  const unsigned i4 , const unsigned i5 , const unsigned i6 , const unsigned i7 )
  {
    assert_shape_bounds( src.shape() , i0 , i1 , i2 , i3 , i4 , i5 , i6 , i7 );

    ViewAssignment< DstViewType >::decrement( dst.m_ptr_on_device );

    dst.m_ptr_on_device =
      src.m_ptr_on_device +
      i7 + src.m_shape.N7 * (
      i6 + src.m_shape.N6 * (
      i5 + src.m_shape.N5 * (
      i4 + src.m_shape.N4 * (
      i3 + src.m_shape.N3 * (
      i2 + src.m_shape.N2 * (
      i1 )))))) + i0 * src.m_stride ;

    ViewAssignment< DstViewType >::increment( dst.m_ptr_on_device );
  }
};

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace KokkosArray */

#endif /* #ifndef KOKKOSARRAY_VIEWRIGHT_HPP */

