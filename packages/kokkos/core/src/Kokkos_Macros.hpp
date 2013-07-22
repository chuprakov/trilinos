/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <KokkosCore_config.h>

#ifndef KOKKOS_MACROS_HPP
#define KOKKOS_MACROS_HPP

namespace Kokkos {
class HostSpace ;
class CudaSpace ;
class Host ;
class Cuda ;
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( __CUDACC__ )

#if ! defined( KOKKOS_HAVE_CUDA )
#error "Compiling Kokkos with Cuda compiler but KOKKOS_ARRAY_CUDA is undefined"
#endif

/*  Compiling with a CUDA compiler for host and device code.
 *
 *  Include <cuda.h> to pick up the CUDA_VERSION macro defined as:
 *    CUDA_VERSION = ( MAJOR_VERSION * 100 ) | ( MINOR_VERSION )
 *
 *  When generating device code the __CUDA_ARCH__ macro is defined as:
 *    __CUDA_ARCH__ = ( MAJOR_CAPABILITY * 100 ) | ( MINOR_CAPABILITY )
 *
 *  Declare functions available on host and device, or device only.
 */

#include <cuda.h>

#if ! defined( CUDA_VERSION )
#error "#include <cuda.h> did not define CUDA_VERSION"
#endif

#if ( CUDA_VERSION < 401 )
#error "Cuda version 4.1 or greater required"
#endif

#if defined( __CUDA_ARCH__ )

#if ( __CUDA_ARCH__ < 200 )
#error "Cuda device capability >= 2.0 is required"
#endif

/*  Compiling with CUDA compiler for device code.
 *  The value of __CUDA_ARCH__ is an integer value XY0
 *  identifying compilation for 'compute_XY' architecture.
 */

namespace Kokkos { typedef CudaSpace ExecutionSpace ; }

#define KOKKOS_FORCEINLINE_FUNCTION  __device__  __host__  __forceinline__
#define KOKKOS_INLINE_FUNCTION  __device__  __host__  inline
#define KOKKOS_FUNCTION         __device__  __host__

#else /* ! #if defined( __CUDA_ARCH__ ) */

/*  Compiling with CUDA compiler for host code. */

namespace Kokkos { typedef HostSpace ExecutionSpace ; }

#define KOKKOS_INLINE_FUNCTION  inline
#define KOKKOS_FUNCTION         /* */

#if defined( __INTEL_COMPILER )
# define KOKKOS_FORCEINLINE_FUNCTION  __forceinline
#else //GCC Compiler
# define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))
#endif

#endif /* ! #if defined( __CUDA_ARCH__ ) */

/* END: #if defined( __CUDACC__ ) */
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#else /* Host memory space and no special markups needed */

namespace Kokkos { typedef HostSpace ExecutionSpace ; }

#define KOKKOS_INLINE_FUNCTION  inline
#define KOKKOS_FUNCTION         /* */

#if defined( __INTEL_COMPILER )
# define KOKKOS_FORCEINLINE_FUNCTION  __forceinline
#else //GCC Compiler
# define KOKKOS_FORCEINLINE_FUNCTION  inline __attribute__((always_inline))
#endif

//----------------------------------------------------------------------------

#if defined( __INTEL_COMPILER )


#if defined( __MIC__ )

/*  Compiling with Intel compiler for execution on an Intel MIC device.
 *  These devices are used in no-offload mode so the Host space is the MIC space.
 */

#endif

#endif

//----------------------------------------------------------------------------

#if defined( __GNUC__ ) /* GNU C   */ || \
    defined( __GNUG__ ) /* GNU C++ */


/*  Compiling with GNU compatible compiler.  */

#endif

//----------------------------------------------------------------------------

#if defined( _OPENMP )

/*  Compiling with in OpenMP mode.
 *  The value of _OPENMP is an integer value YYYYMM
 *  where YYYY and MM are the year and month designation
 *  of the supported OpenMP API version.
 */

#endif /* END: #if defined( _OPENMP ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* END: ! #if defined( __CUDACC__ ) */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#define KOKKOS_RESTRICT_EXECUTION_TO_DATA( DATA_SPACE , DATA_PTR ) \
  Kokkos::VerifyExecutionSpaceCanAccessDataSpace< \
    Kokkos::ExecutionSpace , DATA_SPACE >::verify( DATA_PTR )

#define KOKKOS_RESTRICT_EXECUTION_TO( DATA_SPACE ) \
  Kokkos::VerifyExecutionSpaceCanAccessDataSpace< \
    Kokkos::ExecutionSpace , DATA_SPACE >::verify()




//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_MACROS_HPP */

