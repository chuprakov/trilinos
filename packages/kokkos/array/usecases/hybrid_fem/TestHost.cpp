
/* #define KOKKOS_ARRAY_BOUNDS_CHECK 1 */

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <BoxMeshFixture.hpp>
#include <TestBoxMeshFixture.hpp>
#include <Implicit.hpp>

#include <Kokkos_Host.hpp>

#include <Kokkos_Host_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <TestBoxMeshFixture_macros.hpp>
#include <Implicit_macros.hpp>
#include <SparseLinearSystem_macros.hpp>
#include <SparseLinearSystemFill_macros.hpp>
#include <Kokkos_Clear_macros.hpp>

#include <SparseLinearSystem_Host.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_host_query( comm::Machine machine )
{
  const size_t node_count = Kokkos::Host::detect_node_count();

  std::cout << "P" << comm::rank( machine )
            << ": Host node_count = " << node_count
            << " , node_core_count = "
            << Kokkos::Host::detect_node_core_count();

  std::cout << std::endl ;
}

void test_host_fixture( comm::Machine machine ,
                        size_t numa_node_count ,
                        size_t numa_node_thread_count ,
                        size_t nx , size_t ny , size_t nz )
{
  Kokkos::Host::initialize( numa_node_count , numa_node_thread_count );
  test_box_fixture<Kokkos::Host>( machine , nx , ny , nz );
  Kokkos::Host::finalize();
}

void test_host_implicit( comm::Machine machine , 
                         size_t numa_node_count ,
                         size_t numa_node_thread_count ,
                         size_t node_count_begin ,
                         size_t node_count_end ,
                         size_t count_run )
{
  Kokkos::Host::initialize( numa_node_count , numa_node_thread_count );
  HybridFEM::Implicit::driver<double,Kokkos::Host>( "Host" , machine , node_count_begin , node_count_end , count_run );
  Kokkos::Host::finalize();
}


