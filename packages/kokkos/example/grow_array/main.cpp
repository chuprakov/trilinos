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

#include <iostream>
#include <sstream>

#include <KokkosCore_config.h>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Serial.hpp>
#include <Kokkos_Threads.hpp>
#include <Kokkos_Cuda.hpp>
#include <Kokkos_OpenMP.hpp>

#include <grow_array.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_CUDA )
void grow_array_cuda( int length_array , int span_values );
#endif

int main( int argc , char ** argv )
{
  int num_threads = 4 ;
  int use_numa = 1 ;
  int use_core = 1 ;
  int length_array  = 1000000 ;
  int span_values = 100000000 ;
  

  if ( Kokkos::hwloc::available() ) {
    use_numa = Kokkos::hwloc::get_available_numa_count();
    use_core = Kokkos::hwloc::get_available_cores_per_numa() - 1 ;
    num_threads = use_numa * use_core * Kokkos::hwloc::get_available_threads_per_core();
  }

  {
    std::cout << "Kokkos::Serial" << std::endl ;
    Example::GrowArrayFunctor< Kokkos::Serial >( length_array , span_values );
  }

#if defined( KOKKOS_HAVE_PTHREAD )
  {
    std::cout << "Kokkos::Threads" << std::endl ;
    Kokkos::Threads::initialize( num_threads , use_numa , use_core );
    Example::GrowArrayFunctor< Kokkos::Threads >( length_array , span_values );
    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
  {
    std::cout << "Kokkos::OpenMP" << std::endl ;
    Kokkos::OpenMP::initialize( num_threads , use_numa , use_core );
    Example::GrowArrayFunctor< Kokkos::OpenMP >( length_array , span_values );
    Kokkos::OpenMP::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_CUDA )
  {
    std::cout << "Kokkos::Cuda" << std::endl ;
    grow_array_cuda( length_array , span_values );
  }
#endif

  return 0 ;
}

