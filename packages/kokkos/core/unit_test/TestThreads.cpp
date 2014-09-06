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

#include <Kokkos_Macros.hpp>

#if defined( KOKKOS_HAVE_PTHREAD )

#include <Kokkos_Core.hpp>
#include <Kokkos_hwloc.hpp>

#include <Kokkos_View.hpp>

#include <Kokkos_CrsArray.hpp>

//----------------------------------------------------------------------------

#include <TestViewImpl.hpp>

#include <TestMemoryTracking.hpp>
#include <TestViewAPI.hpp>
#include <TestAggregate.hpp>
#include <TestAtomic.hpp>

#include <TestCrsArray.hpp>
#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestCompilerMacros.hpp>
#include <TestCXX11.hpp>

namespace Test {

class threads : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    // Finalize without initialize is a no-op:
    Kokkos::Threads::finalize();

    const unsigned numa_count       = Kokkos::hwloc::get_available_numa_count();
    const unsigned cores_per_numa   = Kokkos::hwloc::get_available_cores_per_numa();
    const unsigned threads_per_core = Kokkos::hwloc::get_available_threads_per_core();

    unsigned team_count = 0 ;
    unsigned threads_per_team = 0 ;

    // Initialize and finalize with no threads:
    Kokkos::Threads::initialize( 1u );
    Kokkos::Threads::finalize();

    team_count       = std::max( 1u , numa_count );
    threads_per_team = std::max( 2u , cores_per_numa * threads_per_core );

    Kokkos::Threads::initialize( team_count * threads_per_team );
    Kokkos::Threads::finalize();

    team_count       = std::max( 1u , numa_count * 2 );
    threads_per_team = std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );
    Kokkos::Threads::initialize( team_count * threads_per_team );
    Kokkos::Threads::finalize();

    // Quick attempt to verify thread start/terminate don't have race condition:
    team_count       = std::max( 1u , numa_count );
    threads_per_team = std::max( 2u , ( cores_per_numa * threads_per_core ) / 2 );
    for ( unsigned i = 0 ; i < 10 ; ++i ) {
      Kokkos::Threads::initialize( team_count * threads_per_team );
      Kokkos::Threads::sleep();
      Kokkos::Threads::wake();
      Kokkos::Threads::finalize();
    }

    Kokkos::Threads::initialize( team_count * threads_per_team );
    Kokkos::Threads::print_configuration( std::cout );
  }

  static void TearDownTestCase()
  {
    Kokkos::Threads::finalize();
  }
};

TEST_F( threads , init ) {
  ;
}

TEST_F( threads, view_impl) {
  test_view_impl< Kokkos::Threads >();
}

TEST_F( threads, view_api) {
  TestViewAPI< double , Kokkos::Threads >();
}

TEST_F( threads, view_aggregate ) {
  TestViewAggregate< Kokkos::Threads >();
}

TEST_F( threads, long_reduce) {
  TestReduce< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, double_reduce) {
  TestReduce< double ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, team_long_reduce) {
  TestReduceTeam< long ,   Kokkos::Threads >( 100000 );
}

TEST_F( threads, team_double_reduce) {
  TestReduceTeam< double ,   Kokkos::Threads >( 100000 );
}

TEST_F( threads, long_reduce_dynamic ) {
  TestReduceDynamic< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, double_reduce_dynamic ) {
  TestReduceDynamic< double ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, long_reduce_dynamic_view ) {
  TestReduceDynamicView< long ,   Kokkos::Threads >( 1000000 );
}

TEST_F( threads, team_shared_request) {
  TestSharedTeam< Kokkos::Threads >();
}

TEST_F( threads , view_remap )
{
  enum { N0 = 3 , N1 = 2 , N2 = 8 , N3 = 9 };

  typedef Kokkos::View< double*[N1][N2][N3] ,
                             Kokkos::LayoutRight ,
                             Kokkos::Threads > output_type ;

  typedef Kokkos::View< int**[N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Threads > input_type ;

  typedef Kokkos::View< int*[N0][N2][N3] ,
                             Kokkos::LayoutLeft ,
                             Kokkos::Threads > diff_type ;

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

TEST_F( threads , atomics )
{
  const int loop_count = 1e6 ;

  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<unsigned long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<long long int,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<double,Kokkos::Threads>(loop_count,3) ) );

  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,1) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,2) ) );
  ASSERT_TRUE( ( TestAtomic::Loop<float,Kokkos::Threads>(100,3) ) );
}

//----------------------------------------------------------------------------

#if 0
TEST_F( threads , scan_small )
{
  typedef TestScan< Kokkos::Threads , Kokkos::Impl::ThreadsExecUseScanSmall > TestScanFunctor ;
  for ( int i = 0 ; i < 1000 ; ++i ) {
    TestScanFunctor( 10 );
    TestScanFunctor( 10000 );
  }
  TestScanFunctor( 1000000 );
  TestScanFunctor( 10000000 );

  Kokkos::Threads::fence();
}
#endif

TEST_F( threads , scan )
{
  TestScan< Kokkos::Threads >::test_range( 1 , 1000 );
  TestScan< Kokkos::Threads >( 1000000 );
  TestScan< Kokkos::Threads >( 10000000 );
  Kokkos::Threads::fence();
}

//----------------------------------------------------------------------------

TEST_F( threads , team_scan )
{
  TestScanTeam< Kokkos::Threads >( 10 );
  TestScanTeam< Kokkos::Threads >( 10000 );
}

//----------------------------------------------------------------------------

TEST_F( threads , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Threads >() ) );
}


//----------------------------------------------------------------------------
#if defined( KOKKOS_HAVE_CXX11 ) && defined( KOKKOS_HAVE_DEFAULT_DEVICE_TYPE_THREADS )
TEST_F( threads , cxx11 )
{
  if ( Kokkos::Impl::is_same< Kokkos::DefaultExecutionSpace , Kokkos::Threads >::value ) {
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(1) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(2) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(3) ) );
    ASSERT_TRUE( ( TestCXX11::Test< Kokkos::Threads >(4) ) );
  }
}
#endif

} // namespace Test

#endif /* #if defined( KOKKOS_HAVE_PTHREAD ) */
