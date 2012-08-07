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

#include <iostream>

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <KokkosArray_View_macros.hpp> without macros defined"

#else

namespace KokkosArray {

template< class DataType , class LayoutType >
class View< DataType , LayoutType , KOKKOS_MACRO_DEVICE >
  : public Impl::ViewOperator< DataType ,
                               typename LayoutType::array_layout ,
                               KOKKOS_MACRO_DEVICE::memory_space >
{
private:
  template< class , class , class > friend class View ;

  typedef Impl::ViewOperator< DataType ,
                              typename LayoutType::array_layout ,
                              KOKKOS_MACRO_DEVICE::memory_space > oper_type ;

public:

  typedef DataType             data_type ;
  typedef LayoutType           layout_type ;
  typedef KOKKOS_MACRO_DEVICE  device_type ;

  typedef View< data_type , layout_type , device_type >  type ;

  typedef View< typename Impl::add_const< data_type >::type ,
                layout_type , device_type > const_type ;

  typedef View< data_type , layout_type , Host >  HostMirror ;

  typedef typename Impl::remove_all_extents<data_type>::type  value_type ;
  typedef typename LayoutType::array_layout                   array_layout ;
  typedef typename device_type::memory_space                  memory_space ;
  typedef typename device_type::size_type                     size_type ;
  typedef typename oper_type::shape_type                      shape_type ;

  /*------------------------------------------------------------------*/

  static const unsigned Rank = oper_type::shape_type::rank ;

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type rank() const { return oper_type::shape_type::rank ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_0() const { return oper_type::m_shape.N0 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_1() const { return oper_type::m_shape.N1 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_2() const { return oper_type::m_shape.N2 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_3() const { return oper_type::m_shape.N3 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_4() const { return oper_type::m_shape.N4 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_5() const { return oper_type::m_shape.N5 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_6() const { return oper_type::m_shape.N6 ; }

  inline KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension_7() const { return oper_type::m_shape.N7 ; }

  template< typename iType >
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type dimension( const iType & r ) const
    {
      return iType(0) == r ? oper_type::m_shape.N0 : (
             iType(1) == r ? oper_type::m_shape.N1 : (
             iType(2) == r ? oper_type::m_shape.N2 : (
             iType(3) == r ? oper_type::m_shape.N3 : (
             iType(4) == r ? oper_type::m_shape.N4 : (
             iType(5) == r ? oper_type::m_shape.N5 : (
             iType(6) == r ? oper_type::m_shape.N6 : (
             iType(7) == r ? oper_type::m_shape.N7 : 0 )))))));
    }

  /*------------------------------------------------------------------*/
  /** \brief Query shape */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  shape_type shape() const { return oper_type::m_shape ; }

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool () const
  { return 0 != oper_type::m_ptr_on_device ; }

  /** \brief  Query if view to same memory */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const View & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /** \brief  Query if view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const View<rhsDataType,rhsLayout,rhsMemory> & rhs ) const
  {
    return oper_type::m_ptr_on_device == rhs.m_ptr_on_device &&
           oper_type::m_shape == rhs.m_shape ;
  }

  /** \brief  Query if not view to same memory */
  template< class rhsDataType , class rhsLayout , class rhsMemory >
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const View<rhsDataType,rhsLayout,rhsMemory> & rhs ) const
  {
    return oper_type::m_ptr_on_device != rhs.m_ptr_on_device ||
           oper_type::m_shape != rhs.m_shape ;
  }

  /*------------------------------------------------------------------*/

private:

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void internal_private_assign( const shape_type & shape , value_type * ptr )
  {
    oper_type::m_shape          = shape ;
    oper_type::m_ptr_on_device  = ptr ;
    memory_space::increment( oper_type::m_ptr_on_device );
  }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  void internal_private_clear()
  {
    memory_space::decrement( oper_type::m_ptr_on_device );
    oper_type::m_ptr_on_device = 0 ;
  }

public:

  /** \brief  Construct a NULL view */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View() {}

  /** \brief  Construct a view of the array */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View( const View & rhs )
    {
      internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    }

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~View()
  {
    internal_private_clear();
  }

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  View & operator = ( const View & rhs )
  {
    internal_private_clear();
    internal_private_assign( rhs.m_shape , rhs.m_ptr_on_device );
    return *this ;
  }

  /*------------------------------------------------------------------*/
  // Construct a compatible view:

  template< class rhsType , class rhsLayout , class rhsMemory >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space , 
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type > data( rhs.shape() );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  // Assign a compatible view:

  template< class rhsType , class rhsMapSpec , class rhsMemory >
  View & operator = ( const View< rhsType , rhsMapSpec , rhsMemory > & rhs )
    {
      typedef View< rhsType , rhsMapSpec , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space , 
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type > data( rhs.shape() );

      internal_private_clear();
      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );

      return *this ;
    }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a subview */

  template< class rhsType , class rhsLayout , class rhsMemory , class ArgType0 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      typedef Impl::SubShape< shape_type , rhs_shape_type > subshape_type ;

      const subshape_type data( rhs.shape() , arg0 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 , arg2 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape() , arg0 , arg1 , arg2 , arg3 , arg4 , arg5 , arg6 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  template< class rhsType , class rhsLayout , class rhsMemory ,
            class ArgType0 , class ArgType1 , class ArgType2 ,
            class ArgType3 , class ArgType4 , class ArgType5 ,
            class ArgType6 , class ArgType7 >
  View( const View< rhsType , rhsLayout , rhsMemory > & rhs ,
        const ArgType0 & arg0 ,
        const ArgType1 & arg1 ,
        const ArgType2 & arg2 ,
        const ArgType3 & arg3 ,
        const ArgType4 & arg4 ,
        const ArgType5 & arg5 ,
        const ArgType6 & arg6 ,
        const ArgType7 & arg7 )
    {
      typedef View< rhsType , rhsLayout , rhsMemory > rhs_view ;
      typedef typename rhs_view::value_type   rhs_value_type ;
      typedef typename rhs_view::shape_type   rhs_shape_type ;
      typedef typename rhs_view::memory_space rhs_memory_space ;

      typedef typename
        Impl::SubviewAssignable< value_type , memory_space ,
                                 rhs_value_type , rhs_memory_space >::type
          ok_assign ;

      // SubShape<*,*> only exists for valid subshapes:

      const Impl::SubShape< shape_type , rhs_shape_type >
         data( rhs.shape(), arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7 );

      internal_private_assign( data.shape , rhs.m_ptr_on_device + data.offset );
    }

  /*------------------------------------------------------------------*/
  /*------------------------------------------------------------------*/
  /* Creation with allocation of memory on the device */

private:

  /** \brief  This device-specialized method may only be called
   *          within a view constructor.
   */
  void internal_private_create( const std::string & label ,
                                const shape_type shape );

public:

  View( const std::string & label , const shape_type shape )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape );
  }

  explicit View( const std::string & label )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>() );
  }

  View( const std::string & label ,
        const size_t n0 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6) );
  }

  View( const std::string & label ,
        const size_t n0 ,
        const size_t n1 ,
        const size_t n2 ,
        const size_t n3 ,
        const size_t n4 ,
        const size_t n5 ,
        const size_t n6 ,
        const size_t n7 )
  {
    View_create_requires_non_const_data_type< value_type >::success();
    internal_private_create( label , shape_type::template create<memory_space>(n0,n1,n2,n3,n4,n5,n6,n7) );
  }

  /*------------------------------------------------------------------*/
  /** \brief  For unit testing shape mapping. */
  explicit View( const typename oper_type::shape_type & shape )
    { oper_type::m_shape = shape ; }
};

//----------------------------------------------------------------------------

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif

