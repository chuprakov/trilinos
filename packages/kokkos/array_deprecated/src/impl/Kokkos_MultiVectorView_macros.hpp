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

#if ! defined(KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION) || \
    ! defined(KOKKOS_MACRO_DEVICE)                  || \
    ! defined(KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION)

#error "Including <Kokkos_MultiVectorView_macros.hpp> without macros defined"

#else

namespace Kokkos {

//----------------------------------------------------------------------------
/** \brief  Plain-old-data value allocated on a compute device.
 */
template< typename ValueType >
class MultiVectorView< ValueType , KOKKOS_MACRO_DEVICE > {
public:
  typedef ValueType                  value_type ;
  typedef KOKKOS_MACRO_DEVICE        device_type ;
  typedef device_type::memory_space  memory_space ;
  typedef device_type::size_type     size_type ;

  typedef MultiVectorView< value_type , DeviceHost > HostView ;

private:

  enum { Align = Impl::ArrayAlignment< value_type , memory_space >::value };
  enum { RankLength = 0 };
  enum { RankCount  = 1 };
  enum { RankStride = Align <= 1 ? 0 : 2 };

public:

  enum { Contiguous = Align <= 1 };

  /*------------------------------------------------------------------*/
  /** \brief  Query length of vectors */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type length() const { return m_dim[ RankLength ]; }
  
  /** \brief  Query count of vectors */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  size_type count()  const { return m_dim[ RankCount ]; }
  
  /** \brief  Query if NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  operator bool ()  const { return 0 != m_ptr_on_device ; }
  
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator == ( const MultiVectorView & rhs ) const
  {
    return m_ptr_on_device == rhs.m_ptr_on_device &&
           m_dim[ RankCount ] == rhs.m_dim[ RankCount ];
  }
  
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  bool operator != ( const MultiVectorView & rhs ) const
  {
    return m_ptr_on_device != rhs.m_ptr_on_device ||
           m_dim[ RankCount ] != rhs.m_dim[ RankCount ];
  }

  /** \brief  Because memory is contiguous this is exposed */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  value_type * ptr_on_device() const { return m_memory.ptr_on_device(); }

  
  /*------------------------------------------------------------------*/

#if defined(KOKKOS_MACRO_DEVICE_FUNCTION)

  /** \brief  Query value */
  template< typename iTypeP , typename iTypeV >
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP , const iTypeV & iV ) const
    { return m_ptr_on_device[ iP + m_dim[ RankStride ] * iV ]; }
  
  template< typename iTypeP >
  KOKKOS_MACRO_DEVICE_FUNCTION
  value_type & operator()( const iTypeP & iP ) const
    { return m_ptr_on_device[ iP ]; }

#endif /* defined(KOKKOS_MACRO_DEVICE_FUNCTION) */

  /*------------------------------------------------------------------*/
  /** \brief  Construct a NULL view */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MultiVectorView()
    : m_memory(), m_ptr_on_device(0)
    {
      m_dim[ RankLength ] = 0 ;
      m_dim[ RankCount ]  = 0 ;
      m_dim[ RankStride ] = 0 ;
    }

  /** \brief  Construct a view of the array */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MultiVectorView( const MultiVectorView & rhs )
    : m_memory()
    , m_ptr_on_device( rhs.m_ptr_on_device )
    {
      m_dim[ RankLength ] = rhs.m_dim[ RankLength ];
      m_dim[ RankCount ]  = rhs.m_dim[ RankCount ];
      m_dim[ RankStride ] = rhs.m_dim[ RankStride ];
      memory_space::assign_memory_view( m_memory , rhs.m_memory );
    }

  /** \brief  Assign to a view of the rhs.
   *          If the old view is the last view
   *          then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MultiVectorView & operator = ( const MultiVectorView & rhs )
    {
      memory_space::assign_memory_view( m_memory , rhs.m_memory );
      m_ptr_on_device     = rhs.m_ptr_on_device ;
      m_dim[ RankLength ] = rhs.m_dim[ RankLength ];
      m_dim[ RankCount ]  = rhs.m_dim[ RankCount ];
      m_dim[ RankStride ] = rhs.m_dim[ RankStride ];
      return *this ;
    }
  
  /**  \brief  Destroy this view of the value.
   *           If the last view then allocated memory is deallocated.
   */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  ~MultiVectorView()
    {
      memory_space::clear_memory_view( m_memory );
      m_ptr_on_device = 0 ;
      m_dim[ RankLength ] = 0 ;
      m_dim[ RankCount ]  = 0 ;
      m_dim[ RankStride ] = 0 ;
    }

  /*------------------------------------------------------------------*/
  /* \brief  Construct a view to a range of vectors */
  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MultiVectorView( const MultiVectorView & rhs , size_type iBeg ,
                                                 size_type iEnd )
    : m_memory()
    , m_ptr_on_device( 0 )
    {
      if ( iBeg < iEnd && iEnd < rhs.m_dim[ RankCount ] ) {
        m_dim[ RankLength ] = rhs.m_dim[ RankLength ];
        m_dim[ RankCount ]  = iEnd - iBeg ;
        m_dim[ RankStride ] = rhs.m_dim[ RankStride ];
        m_ptr_on_device     = rhs.m_ptr_on_device + m_dim[ RankStride ] * iBeg ;
        memory_space::assign_memory_view( m_memory , rhs.m_memory );
      }
      else {
        m_dim[ RankLength ] = 0 ;
        m_dim[ RankCount ]  = 0 ;
        m_dim[ RankStride ] = 0 ;
        KOKKOS_MACRO_DEVICE_CAN_THROW( Impl::multivector_require_range( iBeg , iEnd , rhs.m_dim[ RankCount] ) );
      }
    }

  inline
  KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION
  MultiVectorView( const MultiVectorView & rhs , size_type iBeg )
    : m_memory()
    , m_ptr_on_device( 0 )
    {
      if ( iBeg < rhs.m_dim[ RankCount ] ) {
        m_dim[ RankLength ] = rhs.m_dim[ RankLength ];
        m_dim[ RankCount ]  = 1 ;
        m_dim[ RankStride ] = rhs.m_dim[ RankStride ];
        m_ptr_on_device     = rhs.m_ptr_on_device + m_dim[ RankStride ] * iBeg ;
        memory_space::assign_memory_view( m_memory , rhs.m_memory );
      }
      else {
        m_dim[ RankLength ] = 0 ;
        m_dim[ RankCount ]  = 0 ;
        m_dim[ RankStride ] = 0 ;
        KOKKOS_MACRO_DEVICE_CAN_THROW( Impl::multivector_require_range( iBeg , iBeg + 1 , rhs.m_dim[ RankCount] ) );
      }
    }

private:

  enum { N = 2 + ( Contiguous ? 0 : 1 ) };

  MemoryView< value_type , memory_space > m_memory ;
  ValueType * m_ptr_on_device ;
  size_type   m_dim[ N ];

  inline
  MultiVectorView( const std::string & label ,
                   size_type arg_length , size_type arg_count )
    : m_memory()
    , m_ptr_on_device( 0 )
    {
      m_dim[ RankLength ] = arg_length ;
      m_dim[ RankCount ]  = arg_count ;

      if ( ! Contiguous ) {
        m_dim[ RankStride ] = m_dim[ RankLength ];
        if ( m_dim[ RankStride ] % Align ) {
          m_dim[ RankStride ] += Align - m_dim[ RankStride ] % Align ;
        }
      }

      memory_space::allocate_memory_view( m_memory , m_dim[ RankStride ] * m_dim[ RankCount ] , label );

      m_ptr_on_device = m_memory.ptr_on_device();
    }

  template< typename V , class M >
  friend
  MultiVectorView< V , M >
  create_labeled_multivector( const std::string & label ,
                              size_t length , size_t count );

  template < typename V , class DeviceDst , class DeviceSrc ,
             bool same_memory_space , bool both_contiguous >
  friend
  class Impl::MultiVectorDeepCopy ;
};

} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* template specialization macros defined */


