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

#include <gtest/gtest.h>

#include <Kokkos_DeviceTPI.hpp>
#include <Kokkos_DeviceTPI_ValueView.hpp>
#include <Kokkos_DeviceTPI_MultiVectorView.hpp>
#include <Kokkos_DeviceTPI_MDArrayView.hpp>
#include <Kokkos_DeviceTPI_ParallelFor.hpp>
#include <Kokkos_DeviceTPI_ParallelReduce.hpp>

//----------------------------------------------------------------------------

#include <Kokkos_DeviceTPI_macros.hpp>

#include <UnitTestDeviceMemoryManagement.hpp>
#include <UnitTestValueView.hpp>
#include <UnitTestMultiVectorView.hpp>
#include <UnitTestMDArrayView.hpp>
#include <UnitTestMDArrayDeepCopy.hpp>
#include <UnitTestMDArrayIndexMap.hpp>
#include <UnitTestReduce.hpp>
#include <UnitTestMultiReduce.hpp>

#include <Kokkos_DeviceClear_macros.hpp>

namespace Test {

class tpi : public ::testing::Test {
  protected:
    static void SetUpTestCase() {
      Kokkos::DeviceTPI::initialize( 4 );
    }
    static void TearDownTestCase() {
      Kokkos::DeviceTPI::finalize();
    }
};



TEST_F( tpi, memory_management_double) {
  UnitTestDeviceMemoryManagement< double, Kokkos::DeviceTPI >();
}

TEST_F( tpi, memory_management_int) {
  UnitTestDeviceMemoryManagement< int, Kokkos::DeviceTPI >();
}

TEST_F( tpi, value_view_double) {
  UnitTestValueView< double, Kokkos::DeviceTPI >();
}

TEST_F( tpi, value_view_int) {
  UnitTestValueView< int, Kokkos::DeviceTPI >();
}

TEST_F( tpi, multi_vector_view_double) {
  UnitTestMultiVectorView< double, Kokkos::DeviceTPI >();
}

TEST_F( tpi, multi_vector_view_int) {
  UnitTestMultiVectorView< int, Kokkos::DeviceTPI >();
}

TEST_F( tpi, mdarray_view_double) {
  UnitTestMDArrayView< double, Kokkos::DeviceTPI >();
}

TEST_F( tpi, mdarray_view_int) {
  UnitTestMDArrayView< int, Kokkos::DeviceTPI >();
}

TEST_F( tpi, mdarray_deep_copy_double) {
  UnitTestMDArrayDeepCopy< double, Kokkos::DeviceTPI >();
}

TEST_F( tpi, mdarray_deep_copy_int) {
  UnitTestMDArrayDeepCopy< int, Kokkos::DeviceTPI >();
}

TEST_F( tpi, mdarray_index_map) {
  UnitTestMDArrayIndexMap< Kokkos::DeviceTPI >();
}

TEST_F( tpi, long_reduce) {
  UnitTestReduce< long ,   Kokkos::DeviceTPI >( 1000000 );
}

TEST_F( tpi, double_reduce) {
  UnitTestReduce< double ,   Kokkos::DeviceTPI >( 1000000 );
}

TEST_F( tpi, long_multi_reduce) {
  UnitTestReduceMulti< long , Kokkos::DeviceTPI >( 1000000 , 7 );
}

} // namespace test

