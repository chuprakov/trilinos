/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
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

#ifndef TEST_AGGREGATE_HPP
#define TEST_AGGREGATE_HPP

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

/*--------------------------------------------------------------------------*/

namespace Test {

struct EmbedArray {};

struct ArrayProxyContiguous {};
struct ArrayProxyStrided {};

template< typename T , unsigned N = 0 , class Proxy = void >
struct Array ;

template< typename T >
struct Array<T,0,ArrayProxyContiguous>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * const value ;
  const unsigned count ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned n ) : value(v), count(n) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,0,Proxy> & rhs ) { return *this ; }
};

template< typename T , unsigned N >
struct Array<T,N,ArrayProxyContiguous>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T * const value ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned ) : value(v) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & rhs ) { return *this ; }
};

template< typename T , unsigned N >
struct Array<T,N,ArrayProxyStrided>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T * const value ;
  const unsigned stride ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned , unsigned s ) : value(v), stride(s) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & rhs ) { return *this ; }
};

template< typename T >
struct Array<T,0,ArrayProxyStrided>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * const value ;
  const unsigned count ;
  const unsigned stride ;

  KOKKOS_INLINE_FUNCTION
  Array( T * v , unsigned n , unsigned s ) : value(v), count(n), stride(s) {}

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,0,Proxy> & rhs ) { return *this ; }
};

template< typename T >
struct Array<T,0,void>
{
public:
  typedef T value_type ;

  enum { StaticLength = 0 };
  T * value ;
  const unsigned count ;

  KOKKOS_INLINE_FUNCTION
  Array() : value(0) , count(0) {}

  template< unsigned N , class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array( const Array<T,N,Proxy> & rhs ) : value(rhs.value), count(N) {}
};

template< typename T , unsigned N >
struct Array<T,N,void>
{
public:
  typedef T value_type ;

  enum { StaticLength = N };
  T value[N] ;

  template< class Proxy >
  KOKKOS_INLINE_FUNCTION
  Array & operator = ( const Array<T,N,Proxy> & ) { return *this ; }
};

} // namespace Test

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template< typename T , unsigned N >
struct AnalyzeShape< Test::Array< T , N > >
  : public ShapeInsert< typename AnalyzeShape< T >::shape , N >::type
{
private:
  typedef AnalyzeShape< T > nested ;
public:

  typedef Test::EmbedArray specialize ;

  typedef typename ShapeInsert< typename nested::shape , N >::type shape ;

  typedef typename nested::array_type   array_type[ N ];
  typedef Test::Array< T , N >          value_type ;
  typedef Test::Array< T , N >          type ;

  typedef const array_type  const_array_type ;
  typedef const value_type  const_value_type ;
  typedef const type        const_type ;

  typedef typename nested::non_const_array_type                    non_const_array_type[ N ];
  typedef Test::Array< typename nested::non_const_value_type , N > non_const_value_type ;
  typedef Test::Array< typename nested::non_const_value_type , N > non_const_type ;
};

template< typename T >
struct AnalyzeShape< Test::Array< T , 0 > >
  : public ShapeInsert< typename AnalyzeShape< T >::shape , 0 >::type
{
private:
  typedef AnalyzeShape< T > nested ;
public:

  typedef Test::EmbedArray specialize ;

  typedef typename ShapeInsert< typename nested::shape , 0 >::type shape ;

  typedef typename nested::array_type * array_type ;
  typedef Test::Array< T , 0 >          value_type ;
  typedef Test::Array< T , 0 >          type ;

  typedef const array_type  const_array_type ;
  typedef const value_type  const_value_type ;
  typedef const type        const_type ;

  typedef typename nested::non_const_array_type  * non_const_array_type ;
  typedef Test::Array< typename nested::non_const_value_type , 0 > non_const_value_type ;
  typedef Test::Array< typename nested::non_const_value_type , 0 > non_const_type ;
};

/*--------------------------------------------------------------------------*/

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ValueType
                     , Test::EmbedArray
                     , LayoutLeft
                     , MemorySpace
                     , MemoryTraits >
{ typedef Test::EmbedArray type ; };

template< class ValueType , class MemorySpace , class MemoryTraits >
struct ViewSpecialize< ValueType
                     , Test::EmbedArray
                     , LayoutRight
                     , MemorySpace
                     , MemoryTraits >
{ typedef Test::EmbedArray type ; };

/*--------------------------------------------------------------------------*/

template<>
struct ViewAssignment< Test::EmbedArray , Test::EmbedArray , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class DT , class DL , class DD , class DM ,
            class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment(       View<DT,DL,DD,DM,Test::EmbedArray> & dst
                , const View<ST,SL,SD,SM,Test::EmbedArray> & src
                , const typename enable_if<(
                    ViewAssignable< ViewTraits<DT,DL,DD,DM> ,
                                    ViewTraits<ST,SL,SD,SM> >::value
                    )>::type * = 0
                  )
  {
    typedef ViewTraits<DT,DL,DD,DM> dst_traits ;
    typedef typename View<DT,DL,DD,DM,ViewDefault>::shape_type   shape_type ;
    typedef typename View<DT,DL,DD,DM,ViewDefault>::stride_type  stride_type ;

    ViewTracking< dst_traits >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_traits >::increment( dst.m_ptr_on_device );
  }
};

template<>
struct ViewAssignment< ViewDefault , Test::EmbedArray , void >
{
  //------------------------------------
  /** \brief  Compatible value and shape */

  template< class ST , class SL , class SD , class SM >
  KOKKOS_INLINE_FUNCTION
  ViewAssignment( typename View<ST,SL,SD,SM,Test::EmbedArray>::array_type & dst
                , const View<ST,SL,SD,SM,Test::EmbedArray> & src
                )
  {
    typedef typename View<ST,SL,SD,SM,Test::EmbedArray>::array_type dst_type ;
    typedef typename dst_type::shape_type   shape_type ;
    typedef typename dst_type::stride_type  stride_type ;

    ViewTracking< dst_type >::decrement( dst.m_ptr_on_device );

    shape_type::assign( dst.m_shape,
                        src.m_shape.N0 , src.m_shape.N1 , src.m_shape.N2 , src.m_shape.N3 ,
                        src.m_shape.N4 , src.m_shape.N5 , src.m_shape.N6 , src.m_shape.N7 );

    stride_type::assign( dst.m_stride , src.m_stride.value );

    dst.m_ptr_on_device = src.m_ptr_on_device ;

    Impl::ViewTracking< dst_type >::increment( dst.m_ptr_on_device );
  }
};


} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

template< class DataType ,
          class Arg1Type ,
          class Arg2Type ,
          class Arg3Type >
class View< DataType , Arg1Type , Arg2Type , Arg3Type , Test::EmbedArray >
  : public ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type >
{
public:

  typedef ViewTraits< DataType , Arg1Type , Arg2Type, Arg3Type > traits ;

private:

  // Assignment of compatible views requirement:
  template< class , class , class , class , class > friend class View ;

  // Assignment of compatible subview requirement:
  template< class , class , class > friend struct Impl::ViewAssignment ;

  typedef Impl::LayoutStride< typename traits::shape_type ,
                              typename traits::array_layout > stride_type ;

  // traits::value_type = Test::Array< T , N >

  typename traits::value_type::value_type * m_ptr_on_device ;
  typename traits::shape_type               m_shape ;
  stride_type                               m_stride ;

public:

  typedef View< typename traits::array_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > array_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > non_const_type ;

  typedef View< typename traits::const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type ,
                typename traits::memory_traits > const_type ;

  typedef View< typename traits::non_const_data_type ,
                typename traits::array_layout ,
                typename traits::device_type::host_mirror_device_type ,
                void > HostMirror ;

  //------------------------------------
  // Shape

  enum { Rank = traits::rank - 1 };

  KOKKOS_INLINE_FUNCTION typename traits::shape_type shape() const { return m_shape ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_0() const { return m_shape.N0 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_1() const { return m_shape.N1 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_2() const { return m_shape.N2 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_3() const { return m_shape.N3 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_4() const { return m_shape.N4 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_5() const { return m_shape.N5 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_6() const { return m_shape.N6 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type dimension_7() const { return m_shape.N7 ; }
  KOKKOS_INLINE_FUNCTION typename traits::size_type size() const
  {
    return   m_shape.N0
           * m_shape.N1
           * m_shape.N2
           * m_shape.N3
           * m_shape.N4
           * m_shape.N5
           * m_shape.N6
           * m_shape.N7
           ;
  }

  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type dimension( const iType & i ) const
    { return Impl::dimension( m_shape , i ); }

  //------------------------------------
  // Destructor, constructors, assignment operators:

  KOKKOS_INLINE_FUNCTION
  ~View() { Impl::ViewTracking< traits >::decrement( m_ptr_on_device ); }

  KOKKOS_INLINE_FUNCTION
  View() : m_ptr_on_device(0)
    {
      traits::shape_type::assign(m_shape,0,0,0,0,0,0,0,0);
      stride_type::assign(m_stride,0);
    }

  KOKKOS_INLINE_FUNCTION
  View( const View & rhs ) : m_ptr_on_device(0)
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
    }

  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize ,
        typename traits::specialize >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Construct or assign compatible view:

  template< class RT , class RL , class RD , class RM , class RS >
  KOKKOS_INLINE_FUNCTION
  View( const View<RT,RL,RD,RM,RS> & rhs )
    : m_ptr_on_device(0)
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize , RS >( *this , rhs );
    }

  template< class RT , class RL , class RD , class RM , class RS >
  KOKKOS_INLINE_FUNCTION
  View & operator = ( const View<RT,RL,RD,RM,RS> & rhs )
    {
      (void) Impl::ViewAssignment<
        typename traits::specialize , RS >( *this , rhs );
      return *this ;
    }

  //------------------------------------
  // Allocation of a managed view with possible alignment padding.

  typedef Impl::if_c< traits::is_managed ,
                      std::string ,
                      Impl::ViewError::allocation_constructor_requires_managed >
   if_allocation_constructor ;

  explicit inline
  View( const typename if_allocation_constructor::type & label ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    {
      //typedef typename traits::device_type   device_type ; // unused
      typedef typename traits::memory_space  memory_space ;
      typedef typename traits::shape_type    shape_type ;
      typedef typename traits::value_type::value_type   scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_shape );

      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );

      Impl::ViewFill< View > init( *this , typename traits::value_type() );
    }

  explicit inline
  View( const AllocateWithoutInitializing & ,
        const typename if_allocation_constructor::type & label ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    {
      //typedef typename traits::device_type   device_type ; // unused
      typedef typename traits::memory_space  memory_space ;
      typedef typename traits::shape_type    shape_type ;
      typedef typename traits::value_type::value_type   scalar_type ;

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_with_padding( m_stride , m_shape );

      m_ptr_on_device = (scalar_type *)
        memory_space::allocate( if_allocation_constructor::select( label ) ,
                                typeid(scalar_type) ,
                                sizeof(scalar_type) ,
                                Impl::capacity( m_shape , m_stride ) );
    }

  //------------------------------------
  // Assign an unmanaged View from pointer, can be called in functors.
  // No alignment padding is performed.

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::value_type::value_type * ,
                      Impl::ViewError::user_pointer_constructor_requires_unmanaged >
    if_user_pointer_constructor ;

  View( typename if_user_pointer_constructor::type ptr ,
        const size_t n0 = 0 ,
        const size_t n1 = 0 ,
        const size_t n2 = 0 ,
        const size_t n3 = 0 ,
        const size_t n4 = 0 ,
        const size_t n5 = 0 ,
        const size_t n6 = 0 ,
        const size_t n7 = 0 )
    : m_ptr_on_device(0)
    {
      typedef typename traits::shape_type   shape_type ;
      //typedef typename traits::value_type::value_type   scalar_type ; // unused

      shape_type ::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_no_padding( m_stride , m_shape );

      m_ptr_on_device = if_user_pointer_constructor::select( ptr );
    }

  //------------------------------------
  // Assign unmanaged View to portion of Device shared memory

  typedef Impl::if_c< ! traits::is_managed ,
                      typename traits::device_type ,
                      Impl::ViewError::device_shmem_constructor_requires_unmanaged >
      if_device_shmem_constructor ;

  explicit KOKKOS_INLINE_FUNCTION
  View( typename if_device_shmem_constructor::type & dev ,
        const unsigned n0 = 0 ,
        const unsigned n1 = 0 ,
        const unsigned n2 = 0 ,
        const unsigned n3 = 0 ,
        const unsigned n4 = 0 ,
        const unsigned n5 = 0 ,
        const unsigned n6 = 0 ,
        const unsigned n7 = 0 )
    : m_ptr_on_device(0)
    {
      typedef typename traits::shape_type   shape_type ;
      typedef typename traits::value_type::value_type   scalar_type ;

      enum { align = 8 };
      enum { mask  = align - 1 };

      shape_type::assign( m_shape, n0, n1, n2, n3, n4, n5, n6, n7 );
      stride_type::assign_no_padding( m_stride , m_shape );

      typedef Impl::if_c< ! traits::is_managed ,
                          scalar_type * ,
                          Impl::ViewError::device_shmem_constructor_requires_unmanaged >
        if_device_shmem_pointer ;

      // Select the first argument:
      m_ptr_on_device = if_device_shmem_pointer::select(
       (scalar_type *) dev.get_shmem( unsigned( sizeof(scalar_type) * Impl::capacity( m_shape , m_stride ) + unsigned(mask) ) & ~unsigned(mask) ) );
    }

  static inline
  unsigned shmem_size( const unsigned n0 = 0 ,
                       const unsigned n1 = 0 ,
                       const unsigned n2 = 0 ,
                       const unsigned n3 = 0 ,
                       const unsigned n4 = 0 ,
                       const unsigned n5 = 0 ,
                       const unsigned n6 = 0 ,
                       const unsigned n7 = 0 )
  {
    enum { align = 8 };
    enum { mask  = align - 1 };

    typedef typename traits::shape_type   shape_type ;
    typedef typename traits::value_type::value_type   scalar_type ;

    shape_type  shape ;
    stride_type stride ;

    traits::shape_type::assign( shape, n0, n1, n2, n3, n4, n5, n6, n7 );
    stride_type::assign_no_padding( stride , shape );

    return unsigned( sizeof(scalar_type) * Impl::capacity( shape , stride ) + unsigned(mask) ) & ~unsigned(mask) ;
  }

  //------------------------------------
  // Is not allocated

  KOKKOS_INLINE_FUNCTION
  bool is_null() const { return 0 == m_ptr_on_device ; }

  //------------------------------------
  // LayoutLeft, rank 2:

  typedef Test::Array< typename traits::value_type::value_type ,
                       traits::value_type::StaticLength ,
                       Test::ArrayProxyStrided > LeftValue ;

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_shape.N1 , m_stride.value );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_shape.N1 , m_stride.value );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< LeftValue , traits, LayoutLeft, 2, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return LeftValue( m_ptr_on_device + i0 , m_shape.N1 , m_stride.value );
    }

  //------------------------------------
  // LayoutRight, rank 2:

  typedef Test::Array< typename traits::value_type::value_type ,
                       traits::value_type::StaticLength ,
                       Test::ArrayProxyContiguous > RightValue ;

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    operator[] ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_shape.N1 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    operator() ( const iType0 & i0 ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_shape.N1 );
    }

  template< typename iType0 >
  KOKKOS_INLINE_FUNCTION
  typename Impl::ViewEnableArrayOper< RightValue , traits, LayoutRight, 2, iType0 >::type
    at( const iType0 & i0 , const int , const int , const int ,
        const int , const int , const int , const int ) const
    {
      KOKKOS_ASSERT_SHAPE_BOUNDS_2( m_shape, i0, 0 );
      KOKKOS_RESTRICT_EXECUTION_TO_DATA( typename traits::memory_space , m_ptr_on_device );

      return RightValue( m_ptr_on_device + i0 , m_shape.N1 );
    }

  //------------------------------------
  // Access to the underlying contiguous storage of this view specialization.
  // These methods are specific to specialization of a view.

  KOKKOS_INLINE_FUNCTION
  typename traits::value_type::value_type * ptr_on_device() const { return m_ptr_on_device ; }

  // Stride of physical storage, dimensioned to at least Rank
  template< typename iType >
  KOKKOS_INLINE_FUNCTION
  void stride( iType * const s ) const
  { Impl::stride( s , m_shape , m_stride ); }

  // Count of contiguously allocated data members including padding.
  KOKKOS_INLINE_FUNCTION
  typename traits::size_type capacity() const
  { return Impl::capacity( m_shape , m_stride ); }
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Test {

template< class DeviceType >
int TestViewAggregate()
{
  typedef Kokkos::View< Test::Array<double,32> * , DeviceType > a32_type ;
  typedef typename a32_type::array_type a32_base_type ;

  typedef Kokkos::View< Test::Array<double> * , DeviceType > a0_type ;
  typedef typename a0_type::array_type a0_base_type ;

  a32_type      a32("a32",100);
  a32_base_type a32_base ;

  a0_type       a0("a0",100,32);
  a0_base_type  a0_base ;

  a32_base = a32 ;
  a0_base = a0 ;

  return 0 ;
}

}

#endif /* #ifndef TEST_AGGREGATE_HPP */
