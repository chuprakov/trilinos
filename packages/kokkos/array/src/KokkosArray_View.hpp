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

#ifndef KOKKOS_VIEW_HPP
#define KOKKOS_VIEW_HPP

#include <string>
#include <impl/KokkosArray_forward.hpp>
#include <impl/KokkosArray_StaticAssert.hpp>
#include <impl/KokkosArray_ArrayTraits.hpp>
#include <impl/KokkosArray_Shape.hpp>

namespace KokkosArray {

//----------------------------------------------------------------------------

template< class DataType , class LayoutType , class DeviceType = LayoutType >
class View {
public:
  typedef DataType    data_type ;
  typedef LayoutType  layout_type ;
  typedef DeviceType  device_type ;

  typedef View< data_type , layout_type , Host >  HostMirror ;

  typedef typename Impl::remove_all_extents<data_type>::type  value_type ;
  typedef typename LayoutType::array_layout                   array_layout ;
  typedef typename device_type::memory_space                  memory_space ;
  typedef typename device_type::size_type                     size_type ;

private:

  typedef typename
    Impl::DefineShape< array_layout , data_type >::type shape_type ;

public:

  typedef Impl::unsigned_< shape_type::rank > Rank ;
 
  size_type rank() const ;

  template< typename iType >
  size_type dimension( const iType & rank ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Value for rank-1 array */
  template< typename iType0 >
  value_type & operator()( const iType0 & i0 ) const ;

  /** \brief  Value for rank-2 array */
  template< typename iType0 , typename iType1 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ) const ;

  /** \brief  Value for rank-3 array */
  template< typename iType0 , typename iType1 , typename iType2 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ) const ;

  /** \brief  Value for rank-4 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ) const ;

  /** \brief  Value for rank-5 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ) const ;

  /** \brief  Value for rank-6 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ) const ;

  /** \brief  Value for rank-7 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ) const ;

  /** \brief  Value for rank-8 array */
  template< typename iType0 , typename iType1 ,
            typename iType2 , typename iType3 ,
            typename iType4 , typename iType5 ,
            typename iType6 , typename iType7 >
  value_type & operator()( const iType0 & i0 ,
                           const iType1 & i1 ,
                           const iType2 & i2 ,
                           const iType3 & i3 ,
                           const iType4 & i4 ,
                           const iType5 & i5 ,
                           const iType6 & i6 ,
                           const iType7 & i7 ) const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  View();

  /** \brief  Construct a view of the array */
  View( const View & rhs );

  /** \brief  Assign to a view of the rhs array.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  View & operator = ( const View & rhs );

  /**  \brief  Destroy this view of the array.
   *           If the last view then allocated memory is deallocated.
   */
  ~View();

  /*------------------------------------------------------------------*/
  /** \brief  Query if NULL view */
  operator bool () const ;

  /** \brief  Query if view to same memory */
  bool operator == ( const View & ) const ;

  /** \brief  Query if not view to same memory */
  bool operator != ( const View & ) const ;
};

//----------------------------------------------------------------------------

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 );
template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 ,
                                    const size_t n6 );

template< class ViewType >
View< typename ViewType::data_type ,
      typename ViewType::layout_type ,
      typename ViewType::device_type >
create( const std::string & label , const size_t n0 , const size_t n1 ,
                                    const size_t n2 , const size_t n3 ,
                                    const size_t n4 , const size_t n5 ,
                                    const size_t n6 , const size_t n7 );

//----------------------------------------------------------------------------

template< class DataTypeDst , class LayoutDst , class DeviceDst ,
          class DataTypeSrc , class LayoutSrc , class DeviceSrc >
void deep_copy( const View<DataTypeDst,LayoutDst,DeviceDst> & dst ,
                const View<DataTypeSrc,LayoutSrc,DeviceSrc> & src );

template< class ScalarType , class LayoutDst , class DeviceDst ,
                             class LayoutSrc , class DeviceSrc >
void deep_copy( const View<ScalarType[],LayoutDst,DeviceDst> & dst ,
                const View<ScalarType[],LayoutSrc,DeviceSrc> & src ,
                const size_t );

template< class ScalarType , class LayoutDst , class DeviceDst >
void deep_copy( const View<ScalarType,LayoutDst,DeviceDst> & dst ,
                const ScalarType & src );

template< class ScalarType , class LayoutSrc , class DeviceSrc >
void deep_copy( ScalarType & dst ,
                const View<ScalarType,LayoutSrc,DeviceSrc> & src );

//----------------------------------------------------------------------------

template < class >
struct template_class_View_requires_a_device_specific_specialization_which_is_not_found ;


template< class DataType , class LayoutType , class DeviceType >
View<DataType,LayoutType,DeviceType>::View()
{
 template_class_View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
}

template< class DataType , class LayoutType , class DeviceType >
View<DataType,LayoutType,DeviceType>::~View()
{
 template_class_View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
}

template< class DataType , class LayoutType , class DeviceType >
View<DataType,LayoutType,DeviceType>::View( const View<DataType,LayoutType,DeviceType> & rhs )
{
 template_class_View_requires_a_device_specific_specialization_which_is_not_found<device_type>::fail();
}

} // namespace KokkosArray

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace KokkosArray {
namespace Impl {

template< class DataType , class ArrayLayout , class MemorySpace >
class ViewOperator ;

}
}

#include <impl/KokkosArray_View_factory.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_VIEW_HPP */

