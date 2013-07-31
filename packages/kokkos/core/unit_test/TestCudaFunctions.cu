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

#include <gtest/gtest.h>

#include <iostream>

#include <Kokkos_Host.hpp>
#include <Kokkos_Cuda.hpp>

#include <Kokkos_View.hpp>
#include <Kokkos_CrsArray.hpp>


//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestViewAPI.hpp>
#include <TestCrsArray.hpp>

#include <TestReduce.hpp>
#include <TestMultiReduce.hpp>

namespace Test {

__global__
void test_abort()
{
  Kokkos::VerifyExecutionSpaceCanAccessDataSpace<
    Kokkos::CudaSpace ,
    Kokkos::HostSpace >::verify();
}


void test_device_cuda_view_impl()
{
  // test_abort<<<32,32>>>(); // Aborts the kernel with CUDA version 4.1 or greater

  test_view_impl< Kokkos::Cuda >();
}

void test_device_cuda_view_api()
{
  TestViewAPI< double , Kokkos::Cuda >();

#if 0
  Kokkos::View<double, Kokkos::Cuda > x("x");
  Kokkos::View<double[1], Kokkos::Cuda > y("y");
  // *x = 10 ;
  // x() = 10 ;
  // y[0] = 10 ;
  // y(0) = 10 ;
#endif
}

void test_device_cuda_crsarray() {
  TestCrsArray< Kokkos::Cuda >();
}

void test_device_cuda_reduce() {
  TestReduce< long ,   Kokkos::Cuda >( 10000000 );
  TestReduce< double , Kokkos::Cuda >( 1000000 );
}

void test_device_cuda_reduce_dynamic() {
  TestReduceDynamic< long ,   Kokkos::Cuda >( 10000000 );
  TestReduceDynamic< double , Kokkos::Cuda >( 1000000 );
}

void test_device_cuda_reduce_dynamic_view() {
  TestReduceDynamicView< long ,   Kokkos::Cuda >( 10000000 );
  TestReduceDynamicView< double , Kokkos::Cuda >( 1000000 );
}

void test_device_cuda_multi_reduce() {
  TestReduceMulti< long , Kokkos::Cuda >( 1000000 , 7 );
}

void test_device_cuda_atomic()
{
  const int loop_count = 1e3 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Cuda>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Cuda>(100,3) ) );
}

}
