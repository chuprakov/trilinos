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

#ifndef KOKKOSARRAY_ARRAYVIEWOPERRIGHT_HPP
#define KOKKOSARRAY_ARRAYVIEWOPERRIGHT_HPP

#include <KokkosArray_Macros.hpp>
#include <KokkosArray_View.hpp>
#include <impl/KokkosArray_Shape.hpp>

namespace KokkosArray {
namespace Impl {

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 1 >
{
private:

  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape , i0 );
      return m_ptr_on_device[ i0 ];
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator[]( const iType0 & i0 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_1( m_shape , i0 );
      return m_ptr_on_device[ i0 ];
    }

  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator * () const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      return ArrayProxyType( m_ptr_on_device );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 2 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape , i0 , i1 );
      return m_ptr_on_device[ i1 + i0 * m_stride ];
    }

  template< typename iType0 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_2( m_shape , i0 , 0 );
      return ArrayProxyType( m_ptr_on_device + i0 * m_stride );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 3 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, i2 );
      return m_ptr_on_device[ i2 + m_shape.N2 * (
                              i1 ) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_3( m_shape, i0, i1, 0 );
      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N2 * ( i1 ) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 4 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, i3 );

      return m_ptr_on_device[ i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ,
                             const iType2 & i2 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_4( m_shape, i0, i1, i2, 0 );

      return ArrayProxyType( m_ptr_on_device + 
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 5 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, i4 );

      return m_ptr_on_device[ i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ,
                             const iType2 & i2 , const iType3 & i3 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_5( m_shape, i0, i1, i2, i3, 0 );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 ))) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 6 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, i5 );

      return m_ptr_on_device[ i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ,
                             const iType2 & i2 , const iType3 & i3 ,
                             const iType4 & i4 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_6( m_shape, i0, i1, i2, i3, i4, 0 );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )))) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 7 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, i6 );

      return m_ptr_on_device[ i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 ))))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ,
                             const iType2 & i2 , const iType3 & i3 ,
                             const iType4 & i4 , const iType5 & i5 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_7( m_shape, i0, i1, i2, i3, i4, i5, 0 );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N6 * ( i5 +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 ))))) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

template< class MemorySpace , typename ScalarType , unsigned Count , class ShapeType>
class ViewOper< MemorySpace ,
                Array< ScalarType , Count , void > ,
                ShapeType , LayoutRight , 8 >
{
private:
  template< class , class , class , class , class > friend class KokkosArray::View ;
  template< class , class , class > friend class KokkosArray::Impl::ViewAssignment ;

  ScalarType * m_ptr_on_device ;
  ShapeType    m_shape ;
  unsigned     m_stride ;

  typedef Array< ScalarType , Count , ArrayProxyContiguous > ArrayProxyType ;

public:

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  KOKKOSARRAY_INLINE_FUNCTION
  ScalarType & operator()( const iType0 & i0 , const iType1 & i1 ,
                          const iType2 & i2 , const iType3 & i3 ,
                          const iType4 & i4 , const iType5 & i5 ,
                          const iType6 & i6 , const iType7 & i7 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, i7 );

      return m_ptr_on_device[ i7 + m_shape.N7 * (
                              i6 + m_shape.N6 * (
                              i5 + m_shape.N5 * (
                              i4 + m_shape.N4 * (
                              i3 + m_shape.N3 * (
                              i2 + m_shape.N2 * (
                              i1 )))))) + i0 * m_stride ];
    }

  template< typename iType0 , typename iType1 , typename iType2 ,
            typename iType3 , typename iType4 , typename iType5 ,
            typename iType6 >
  KOKKOSARRAY_INLINE_FUNCTION
  ArrayProxyType operator()( const iType0 & i0 , const iType1 & i1 ,
                             const iType2 & i2 , const iType3 & i3 ,
                             const iType4 & i4 , const iType5 & i5 ,
                             const iType6 & i6 ) const
    {
      KOKKOSARRAY_RESTRICT_EXECUTION_TO_DATA( MemorySpace , m_ptr_on_device );
      KOKKOSARRAY_ASSERT_SHAPE_BOUNDS_8( m_shape, i0, i1, i2, i3, i4, i5, i6, 0 );

      return ArrayProxyType( m_ptr_on_device +
                             m_shape.N7 * ( i6 +
                             m_shape.N6 * ( i5 +
                             m_shape.N5 * ( i4 +
                             m_shape.N4 * ( i3 +
                             m_shape.N3 * ( i2 +
                             m_shape.N2 * ( i1 )))))) +
                             m_stride * i0 );
    }

};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOSARRAY_ARRAYVIEWOPERRIGHT_HPP */

