/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Kokkos_ValueView_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType >
class ValueView< ValueType , KOKKOS_MACRO_DEVICE > {
public:
  typedef ValueType                  value_type ;
  typedef KOKKOS_MACRO_DEVICE        device_type ;
  typedef device_type::memory_space  memory_space ;

  typedef ValueView< value_type , DeviceHost >  HostView ;

  /*------------------------------------------------------------------*/

#if defined(KOKKOS_MACRO_DEVICE_FUNCTION)

  /** \brief  Access value */
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator* () const 
  { return * m_memory.ptr_on_device(); }

#endif /* defined(KOKKOS_MACRO_DEVICE_FUNCTION) */

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool() const
  { return m_memory.operator bool(); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const ValueView & rhs ) const
  { return m_memory.operator==( rhs.m_memory ); }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const ValueView & rhs ) const
  { return m_memory.operator!=( rhs.m_memory ); }

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ValueView() : m_memory() {}

  /** \brief  Construct a view of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ValueView( const ValueView & rhs )
    : m_memory()
    { memory_space::assign_memory_view( m_memory , rhs.m_memory ); }

  /** \brief  Assign to a view of the rhs.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ValueView & operator = ( const ValueView & rhs )
    { memory_space::assign_memory_view( m_memory , rhs.m_memory ); return *this ; }
  
  /**  \brief  Destroy this view of the value.
   *           If the last view then allocated memory is deallocated.
   */
  ~ValueView()
    { memory_space::clear_memory_view( m_memory ); }

  /*------------------------------------------------------------------*/

  template< class DeviceRHS >
  inline
  explicit
  ValueView( const ValueView< value_type , DeviceRHS > & rhs )
    : m_memory()
    { memory_space::assign_memory_view( m_memory , rhs.m_memory ); }

  /*------------------------------------------------------------------*/

  inline
  bool operator == ( const ValueView & rhs )
  { return m_memory == rhs.m_memory(); }

  template< class DeviceRHS >
  inline
  bool operator != ( const ValueView & rhs )
  { return m_memory != rhs.m_memory(); }

  /*------------------------------------------------------------------*/

  template< class DeviceRHS >
  inline
  bool operator == ( const ValueView< value_type , DeviceRHS > & rhs )
  { return m_memory == rhs.m_memory(); }

  template< class DeviceRHS >
  inline
  bool operator != ( const ValueView< value_type , DeviceRHS > & rhs )
  { return m_memory != rhs.m_memory(); }

  /*------------------------------------------------------------------*/

#if defined(KOKKOS_MACRO_DEVICE_FUNCTION)

  /** \brief  Allow the ValueView to be a parallel reduce
   *          'finalize functor' that assigns the reduced value
   *          on the device.
   */
  inline
  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( const value_type & rhs ) const
    { * m_memory.ptr_on_device() = rhs ; }

#endif /* defined(KOKKOS_MACRO_DEVICE_FUNCTION) */

private:

  MemoryView< value_type , memory_space > m_memory ;

  inline
  explicit ValueView( const std::string & label )
    : m_memory()
    { memory_space::allocate_memory_view( m_memory , 1 , label ); }

  template< typename V , class D >
  friend
  ValueView< V , D >
  create_labeled_value( const std::string & label );

  template< typename , class , class >
  friend
  class Impl::ValueDeepCopy ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {



} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
#endif /* KOKKOS_VALUEVIEW_HPP */


