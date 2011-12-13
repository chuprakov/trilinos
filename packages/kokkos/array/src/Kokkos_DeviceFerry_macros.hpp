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

#if ! defined( KOKKOS_DEVICE_FERRY ) || \
    defined( KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION ) || \
    defined( KOKKOS_MACRO_DEVICE ) || \
    defined( KOKKOS_MACRO_DEVICE_FUNCTION ) || \
    defined( KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION )

#error "Including <Kokkos_DeviceFerry_macros.hpp> with macros already defined"

#else

/*  If compiling with Ferry compiler
 *  then must clarify what functions are only available on the device
 *  versus available on both the device and host.
 */

#define KOKKOS_MACRO_DEVICE_TEMPLATE_SPECIALIZATION /* */
#define KOKKOS_MACRO_DEVICE                      Kokkos::DeviceFerry
#define KOKKOS_MACRO_DEVICE_MEMORY               Kokkos::DeviceFerry
#define KOKKOS_MACRO_DEVICE_AND_HOST_FUNCTION    __declspec(target(mic))
#define KOKKOS_MACRO_DEVICE_FUNCTION_FERRY 		 _Shared
#define KOKKOS_MACRO_DEVICE_FUNCTION             __declspec(target(mic))

/* Device-only functions are not compiled unless compiling for the device. */
#if defined( __MIC__ )
#define KOKKOS_MACRO_DEVICE_CAN_THROW( expr )  /* */
#define KOKKOS_MACRO_CHECK( expr )             /* */
#else
#define KOKKOS_MACRO_DEVICE_CAN_THROW( expr ) expr
#define KOKKOS_MACRO_CHECK( expr )  /* */
#endif

#endif

