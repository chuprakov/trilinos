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

#ifndef KOKKOS_CUDASPACE_HPP
#define KOKKOS_CUDASPACE_HPP

#if defined( KOKKOS_HAVE_CUDA )

//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

/*  Compiling with a CUDA compiler.
 *
 *  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
 *    CUDA_VERSION = ( MAJOR_VERSION * 1000 ) + ( MINOR_VERSION * 10 )
 *
 *  When generating device code the __CUDA_ARCH__ macro is defined as:
 *    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) + ( MINOR_CAPABILITY * 10 )
 */

#include <cuda_runtime.h>
#include <cuda.h>

#if ! defined( CUDA_VERSION )
#error "#include <cuda.h> did not define CUDA_VERSION"
#endif

#if ( CUDA_VERSION < 4010 )
#error "Cuda version 4.1 or greater required"
#endif

#if defined( __CUDA_ARCH__ ) && ( __CUDA_ARCH__ < 200 )
/*  Compiling with CUDA compiler for device code. */
#error "Cuda device capability >= 2.0 is required"
#endif

#endif /* #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------

#include <iosfwd>
#include <typeinfo>
#include <string>

#include <Kokkos_Macros.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Cuda/Kokkos_Cuda_abort.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda on-device memory management */

class CudaSpace {
public:

  typedef Impl::MemorySpaceTag  kokkos_tag ;
  typedef CudaSpace             memory_space ;
  typedef unsigned int          size_type ;
  typedef Kokkos::Cuda          execution_space ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device
   *          with size = scalar_size * scalar_count.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

  /** \brief  Increment the reference count of the block of memory
   *          in which the input pointer resides.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void increment( const void * );

  /** \brief  Decrement the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count falls to zero the memory is deallocated.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void decrement( const void * );

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /*--------------------------------*/
  /** \brief  Cuda specific function to attached texture object to an allocation.
   *          Output the texture object, base pointer, and offset from the input pointer.
   */
#if defined( __CUDACC__ )
  static void texture_object_attach( const void            * const arg_ptr
                                   , ::cudaChannelFormatDesc const & arg_desc
                                   , ::cudaTextureObject_t * const arg_tex_obj
                                   , void const           ** const arg_alloc_ptr
                                   , int                   * const arg_offset
                                   );
#endif

  /*--------------------------------*/
  /** \brief  Error reporting for HostSpace attempt to access CudaSpace */
  static void access_error();
  static void access_error( const void * const );
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Cuda memory that is accessible to Host execution space
 *          through Cuda's unified virtual memory (UVM) runtime.
 */
class CudaUVMSpace {
public:

  typedef Impl::MemorySpaceTag  kokkos_tag ;
  typedef CudaUVMSpace          memory_space ;
  typedef unsigned int          size_type ;
  typedef Cuda                  execution_space ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device
   *          with size = scalar_size * scalar_count.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

  /** \brief  Increment the reference count of the block of memory
   *          in which the input pointer resides.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void increment( const void * );

  /** \brief  Decrement the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count falls to zero the memory is deallocated.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void decrement( const void * );

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );

  /** \brief  Cuda specific function to attached texture object to an allocation.
   *          Output the texture object, base pointer, and offset from the input pointer.
   */
#if defined( __CUDACC__ )
  static void texture_object_attach( const void            * const arg_ptr
                                   , ::cudaChannelFormatDesc const & arg_desc
                                   , ::cudaTextureObject_t * const arg_tex_obj
                                   , void const           ** const arg_alloc_ptr
                                   , int                   * const arg_offset
                                   );
#endif
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Host memory that is accessible to Cuda execution space
 *          through Cuda's host-pinned memory allocation.
 */
class CudaHostPinnedSpace {
public:

  typedef Impl::MemorySpaceTag        kokkos_tag ;
  typedef CudaHostPinnedSpace         memory_space ;
  typedef unsigned int                size_type ;

  /** \brief  Memory is in HostSpace so use the HostSpace::execution_space */
  typedef HostSpace::execution_space  execution_space ;

  /** \brief  Allocate a contiguous block of memory on the Cuda device
   *          with size = scalar_size * scalar_count.
   *
   *  The input label is associated with the block of memory.
   *  The block of memory is tracked via reference counting where
   *  allocation gives it a reference count of one.
   *
   *  Allocation may only occur on the master thread of the process.
   */
  static void * allocate( const std::string    & label ,
                          const std::type_info & scalar_type ,
                          const size_t           scalar_size ,
                          const size_t           scalar_count );

  /** \brief  Increment the reference count of the block of memory
   *          in which the input pointer resides.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void increment( const void * );

  /** \brief  Decrement the reference count of the block of memory
   *          in which the input pointer resides.  If the reference
   *          count falls to zero the memory is deallocated.
   *
   *          Reference counting only occurs on the master thread.
   */
  static void decrement( const void * );

  /** \brief  Print all tracked memory to the output stream. */
  static void print_memory_view( std::ostream & );

  /** \brief  Retrieve label associated with the input pointer */
  static std::string query_label( const void * );
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template<> struct DeepCopy< CudaSpace , CudaSpace > { DeepCopy( void * dst , const void * src , size_t ); };
template<> struct DeepCopy< CudaSpace , HostSpace > { DeepCopy( void * dst , const void * src , size_t ); };
template<> struct DeepCopy< HostSpace , CudaSpace > { DeepCopy( void * dst , const void * src , size_t ); };

template<> struct DeepCopy< CudaSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< CudaUVMSpace , CudaSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaUVMSpace , HostSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< CudaSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< CudaHostPinnedSpace , CudaSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};

template<> struct DeepCopy< CudaHostPinnedSpace , HostSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};


template<> struct DeepCopy< HostSpace , CudaUVMSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , CudaSpace >( dst , src , n ); }
};

template<> struct DeepCopy< HostSpace , CudaHostPinnedSpace >
{
  inline
  DeepCopy( void * dst , const void * src , size_t n )
  { (void) DeepCopy< HostSpace , HostSpace >( dst , src , n ); }
};

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** Running in CudaSpace attempting to access HostSpace: error */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::HostSpace >
{
  KOKKOS_INLINE_FUNCTION static void verify( void )
    { Kokkos::cuda_abort("Cuda code attempted to access HostSpace memory"); }

  KOKKOS_INLINE_FUNCTION static void verify( const void * )
    { Kokkos::cuda_abort("Cuda code attempted to access HostSpace memory"); }
};

/** Running in CudaSpace accessing CudaUVMSpace: ok */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::CudaUVMSpace >
{
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

/** Running in CudaSpace accessing CudaHostPinnedSpace: ok */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::CudaSpace , Kokkos::CudaHostPinnedSpace >
{
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

/** Running in CudaSpace attempting to access an unknown space: error */
template< class OtherSpace >
struct VerifyExecutionCanAccessMemorySpace<
  typename enable_if< ! is_same<Kokkos::CudaSpace,OtherSpace>::value , Kokkos::CudaSpace >::type ,
  OtherSpace >
{
  KOKKOS_INLINE_FUNCTION static void verify( void )
    { Kokkos::cuda_abort("Cuda code attempted to access unknown Space memory"); }

  KOKKOS_INLINE_FUNCTION static void verify( const void * )
    { Kokkos::cuda_abort("Cuda code attempted to access unknown Space memory"); }
};

//----------------------------------------------------------------------------
/** Running in HostSpace attempting to access CudaSpace */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaSpace >
{
#if defined( KOKKOS_USE_CUDA_UVM )
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
#else
  inline static void verify( void ) { CudaSpace::access_error(); }
  inline static void verify( const void * p ) { CudaSpace::access_error(p); }
#endif
};

/** Running in HostSpace accessing CudaUVMSpace is OK */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaUVMSpace >
{
  inline static void verify( void ) { }
  inline static void verify( const void * ) { }
};

/** Running in HostSpace accessing CudaHostPinnedSpace is OK */
template<>
struct VerifyExecutionCanAccessMemorySpace< Kokkos::HostSpace , Kokkos::CudaHostPinnedSpace >
{
  KOKKOS_INLINE_FUNCTION static void verify( void ) {}
  KOKKOS_INLINE_FUNCTION static void verify( const void * ) {}
};

//----------------------------------------------------------------------------

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_HAVE_CUDA ) */
#endif /* #define KOKKOS_CUDASPACE_HPP */

