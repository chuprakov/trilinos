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

#ifndef KOKKOS_VALUEVIEW_HPP
#define KOKKOS_VALUEVIEW_HPP

#include <cstddef>
#include <string>
#include <Kokkos_ArrayForwardDeclarations.hpp>

namespace Kokkos {

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
ValueView< ValueType , DeviceType >
create_labeled_value( const std::string & label );

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType , class DeviceType >
class ValueView {
public:
  typedef ValueType  value_type ;
  typedef DeviceType device_type ;

  /*------------------------------------------------------------------*/
  /** \brief  Access value */
  value_type & operator* () const ;

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  ValueView();

  /** \brief  Construct a view of the array */
  ValueView( const ValueView & rhs );

  /** \brief  Assign to a view of the rhs.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  ValueView & operator = ( const ValueView & rhs );
  
  /**  \brief  Destroy this view of the value.
   *           If the last view then allocated memory is deallocated.
   */
  ~ValueView();

  /*------------------------------------------------------------------*/
  /** \brief  Allow the ValueView to be a parallel reduce
   *          'finalize functor' that assigns the reduced value
   *          on the device.
   */
  void operator()( const value_type & rhs ) const ;

private:

  explicit ValueView( const std::string & label );

  template< typename V , class D >
  friend
  ValueView< V , D >
  create_labeled_value( const std::string & label );
};

//----------------------------------------------------------------------------

template< typename ValueType , class DeviceType >
inline
ValueView< ValueType , DeviceType >
create_labeled_value( const std::string & label )
{ return ValueView< ValueType , DeviceType >( label ); }

template< typename ValueType , class DeviceType >
inline
ValueView< ValueType , DeviceType >
create_value()
{ return create_labeled_value< ValueType , DeviceType >( std::string() ); }

template< typename ViewType >
inline
ValueView< typename ViewType::value_type , 
           typename ViewType::device_type >
create_labeled_value( const std::string & label )
{
  return create_labeled_value< typename ViewType::value_type , 
                               typename ViewType::device_type >( label );
}

template< typename ViewType >
inline
ValueView< typename ViewType::value_type , 
           typename ViewType::device_type >
create_value( const std::string & label )
{
  return create_labeled_value< typename ViewType::value_type , 
                               typename ViewType::device_type >
    ( std::string() );
}

} // namespace Kokkos

//----------------------------------------------------------------------------
// Partial specializations for known devices

#if defined( KOKKOS_DEVICE_HOST )
#include <impl/Kokkos_DeviceHost_macros.hpp>
#include <impl/Kokkos_ValueView_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_TPI )
#include <impl/Kokkos_DeviceTPI_macros.hpp>
#include <impl/Kokkos_ValueView_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif

#if defined( KOKKOS_DEVICE_CUDA )
#include <impl/Kokkos_DeviceCuda_macros.hpp>
#include <impl/Kokkos_ValueView_macros.hpp>
#include <impl/Kokkos_DeviceClear_macros.hpp>
#endif

//----------------------------------------------------------------------------

#endif /* KOKKOS_VALUEVIEW_HPP */


