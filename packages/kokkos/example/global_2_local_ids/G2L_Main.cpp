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

#include <Kokkos_Serial.hpp>
#include <Kokkos_Threads.hpp>

#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif

#ifdef KOKKOS_HAVE_CUDA
#include <Kokkos_Cuda.hpp>
#endif

#include <Kokkos_hwloc.hpp>


#include <G2L.hpp>

namespace G2L {

size_t run_serial(unsigned num_ids, unsigned num_find_iterations)
{
  std::cout << "Serial" << std::endl;
  return run_test<Kokkos::Serial>(num_ids,num_find_iterations);
}

size_t run_threads(unsigned num_ids, unsigned num_find_iterations)
{
  std::cout << "Threads" << std::endl;
  return run_test<Kokkos::Threads>(num_ids,num_find_iterations);
}

size_t run_openmp(unsigned num_ids, unsigned num_find_iterations)
{
#ifdef KOKKOS_HAVE_OPENMP
  std::cout << "OpenMP" << std::endl;
  return run_test<Kokkos::OpenMP>(num_ids,num_find_iterations);
#else
  return 0;
#endif
}

#ifdef KOKKOS_HAVE_CUDA
extern size_t run_cuda(unsigned num_ids, unsigned num_find_iterations);
#else
size_t run_cuda(unsigned num_ids, unsigned num_find_iterations)
{
  return 0;
}
#endif

} // namespace G2L


int main(int argc, char *argv[])
{
  unsigned num_ids = 100000;
  unsigned num_find_iterations = 1000;

  if (argc == 3) {
    num_ids = atoi(argv[1]);
    num_find_iterations = atoi(argv[2]);
  }
  else if (argc != 1) {
    std::cout << argv[0] << " num_ids num_find_iterations" << std::endl;
    return 0;
  }


  // query the topology of the host
  unsigned team_count = 1 ;
  unsigned threads_count = 4 ;

  //avoid unused variable warning
  (void)team_count;

  if (Kokkos::hwloc::available()) {
    threads_count = Kokkos::hwloc::get_available_numa_count() *
                    Kokkos::hwloc::get_available_cores_per_numa();
  }

  std::cout << "Threads: " << threads_count << std::endl;
  std::cout << "Number of ids: " << num_ids << std::endl;
  std::cout << "Number of find iterations: " << num_find_iterations << std::endl;

  size_t num_errors = 0;

  num_errors += G2L::run_serial(num_ids,num_find_iterations);

#ifdef KOKKOS_HAVE_PTHREAD
  Kokkos::Threads::initialize( threads_count );
  num_errors += G2L::run_threads(num_ids,num_find_iterations);
  Kokkos::Threads::finalize();
#endif

#ifdef KOKKOS_HAVE_OPENMP
  Kokkos::OpenMP::initialize( threads_count );
  num_errors += G2L::run_openmp(num_ids,num_find_iterations);
  Kokkos::OpenMP::finalize();
#endif

#ifdef KOKKOS_HAVE_CUDA
  Kokkos::Cuda::host_mirror_device_type::initialize(1);
  Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
  num_errors += G2L::run_cuda(num_ids,num_find_iterations);
  Kokkos::Cuda::finalize();
#endif

  return num_errors;
}

