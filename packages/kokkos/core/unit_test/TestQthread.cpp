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

#include <Kokkos_Atomic.hpp>
#include <Kokkos_Qthread.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_CrsArray.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestMemoryTracking.hpp>
#include <TestViewAPI.hpp>

#include <TestCrsArray.hpp>
#include <TestTeam.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestMultiReduce.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>

namespace Test {

class qthread : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    int threads_count = 4 ;
    Kokkos::Qthread::initialize( threads_count );
    Kokkos::Qthread::print_configuration( std::cout , true );
  }

  static void TearDownTestCase()
  {
    Kokkos::Qthread::finalize();
  }
};

TEST_F( qthread , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Qthread >() ) );
}

TEST_F( qthread, view_impl) {
  test_view_impl< Kokkos::Qthread >();
}

TEST_F( qthread, view_api) {
  TestViewAPI< double , Kokkos::Qthread >();
}


TEST_F( qthread, crsarray) {
  TestCrsArray< Kokkos::Qthread >();
}

TEST_F( qthread, long_reduce) {
  TestReduce< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, double_reduce) {
  TestReduce< double ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Qthread >( 1000000 );
}

TEST_F( qthread, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Qthread >( 1000000 );
}


TEST_F( qthread , atomics )
{
  const int loop_count = 1e4 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Qthread>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Qthread>(100,3) ) );
}

TEST_F( qthread , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Qthread > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Qthread > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Qthread > diff_type ;

  output_type output( "output" , N0 );
  input_type  input ( "input" , N0 , N1 );
  diff_type   diff  ( "diff" , N0 );

  int value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    input(i0,i1,i2,i3) = ++value ;
  }}}}

  // Kokkos::deep_copy( diff , input ); // throw with incompatible shape
  Kokkos::deep_copy( output , input );

  value = 0 ;
  for ( size_t i3 = 0 ; i3 < N3 ; ++i3 ) {
  for ( size_t i2 = 0 ; i2 < N2 ; ++i2 ) {
  for ( size_t i1 = 0 ; i1 < N1 ; ++i1 ) {
  for ( size_t i0 = 0 ; i0 < N0 ; ++i0 ) {
    ++value ;
    ASSERT_EQ( value , ((int) output(i0,i1,i2,i3) ) );
  }}}}
}

//----------------------------------------------------------------------------


TEST_F( qthread , view_aggregate )
{
  TestViewAggregate< Kokkos::Qthread >();
}

//----------------------------------------------------------------------------

TEST_F( qthread , scan )
{
  for ( int i = 0 ; i < 1000 ; ++i ) {
    TestScan< Kokkos::Qthread >( 10 );
    TestScan< Kokkos::Qthread >( 10000 );
  }
  TestScan< Kokkos::Qthread >( 1000000 );
  TestScan< Kokkos::Qthread >( 10000000 );
  Kokkos::Qthread::fence();
}

TEST_F( qthread, team_shared ) {
  TestSharedTeam< Kokkos::Qthread >();
}

TEST_F( qthread , team_scan )
{
  TestScanTeam< Kokkos::Qthread >( 10 );
//  TestScanTeam< Kokkos::Qthread >( 10000 );
}

//----------------------------------------------------------------------------

} // namespace test

