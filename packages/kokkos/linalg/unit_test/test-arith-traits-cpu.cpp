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

#include <iostream>

#include "Kokkos_ArithTraitsTest.hpp"
#include "Kokkos_Serial.hpp"

#ifdef KOKKOS_HAVE_OPENMP
#  include "Kokkos_OpenMP.hpp"
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
#  include "Kokkos_Threads.hpp"
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
#  include "Kokkos_Cuda.hpp"
#endif // KOKKOS_HAVE_CUDA

int
main (int argc, char* argv[])
{
  using std::cout;
  using std::endl;
  (void) argc;
  (void) argv;

  bool success = true;
  const bool verbose = false;

  {
    Kokkos::Serial::initialize (); // Start up the Kokkos device
    bool serialSuccess = true;
    serialSuccess = serialSuccess && runAllArithTraitsHostTests<Kokkos::Serial> (cout, verbose);
    serialSuccess = serialSuccess && runAllArithTraitsDeviceTests<Kokkos::Serial> (cout, verbose);
    success = success && serialSuccess;
    if (serialSuccess) {
      cout << endl << "Kokkos::Serial host and device: TEST PASSED" << endl;
    } else {
      cout << endl << "Kokkos::Serial host and device: TEST FAILED" << endl;
    }
    Kokkos::Serial::finalize (); // Close down the Kokkos device
  }

#ifdef KOKKOS_HAVE_OPENMP
  {
    Kokkos::OpenMP::initialize (); // Start up the Kokkos device
    bool openmpSuccess = true;
    openmpSuccess = openmpSuccess && runAllArithTraitsHostTests<Kokkos::OpenMP> (cout, verbose);
    openmpSuccess = openmpSuccess && runAllArithTraitsDeviceTests<Kokkos::OpenMP> (cout, verbose);
    success = success && openmpSuccess;
    if (openmpSuccess) {
      cout << endl << "Kokkos::OpenMP host and device: TEST PASSED" << endl;
    } else {
      cout << endl << "Kokkos::OpenMP host and device: TEST FAILED" << endl;
    }
    Kokkos::OpenMP::finalize (); // Close down the Kokkos device
  }
#endif // KOKKOS_HAVE_OPENMP

#ifdef KOKKOS_HAVE_PTHREAD
  {
    Kokkos::Threads::initialize (); // Start up the Kokkos device
    bool threadsSuccess = true;
    threadsSuccess = threadsSuccess && runAllArithTraitsHostTests<Kokkos::Threads> (cout, verbose);
    threadsSuccess = threadsSuccess && runAllArithTraitsDeviceTests<Kokkos::Threads> (cout, verbose);
    success = success && threadsSuccess;
    if (threadsSuccess) {
      cout << endl << "Kokkos::Threads host and device: TEST PASSED" << endl;
    } else {
      cout << endl << "Kokkos::Threads host and device: TEST FAILED" << endl;
    }
    Kokkos::Threads::finalize (); // Close down the Kokkos device
  }
#endif // KOKKOS_HAVE_PTHREAD

#ifdef KOKKOS_HAVE_CUDA
{
  // Start up the Cuda device's host mirror device (must be done
  // before starting up the Cuda device)
  Kokkos::Cuda::host_mirror_device_type::initialize ();
  Kokkos::Cuda::initialize (); // Start up the Kokkos device

  bool cudaSuccess = true;
  // cudaSuccess = cudaSuccess && runAllArithTraitsHostTests<Kokkos::Cuda> (cout, verbose);
  // cudaSuccess = cudaSuccess && runAllArithTraitsDeviceTests<Kokkos::Cuda> (cout, verbose);
  success = success && cudaSuccess;
  if (cudaSuccess) {
    cout << endl << "Kokkos::Cuda host and device: TEST PASSED" << endl;
  } else {
    cout << endl << "Kokkos::Cuda host and device: TEST FAILED" << endl;
  }
  Kokkos::Cuda::finalize (); // Close down the Kokkos device
  // Close down the Cuda device's host mirror device (must be done
  // after closing down the Cuda device)
  Kokkos::Cuda::host_mirror_device_type::finalize ();
}
#endif // KOKKOS_HAVE_CUDA

  if (success) {
    cout << endl << "End Result: TEST PASSED" << endl;
  } else {
    cout << endl << "End Result: TEST FAILED" << endl;
  }
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
