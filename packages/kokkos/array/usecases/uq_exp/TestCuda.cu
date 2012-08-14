
/* #define KOKKOS_ARRAY_BOUNDS_CHECK 1 */

#include <iostream>
#include <stdexcept>
#include <limits>
#include <utility>

#include <BoxMeshFixture.hpp>
#include <Explicit.hpp>

#include <KokkosArray_Cuda.hpp>

#include <KokkosArray_Cuda_macros.hpp>
#include <ParallelDataMap_macros.hpp>
#include <Explicit_macros.hpp>
#include <KokkosArray_Clear_macros.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

void test_cuda_explicit( comm::Machine machine , 
                         size_t elem_count_begin ,
                         size_t elem_count_end ,
                         size_t uq_count_begin ,
                         size_t uq_count_end ,
                         size_t count_run )
{
  const size_t comm_rank = comm::rank( machine );
  const size_t dev_count = KokkosArray::Cuda::detect_device_count();
  const size_t dev_rank  = comm_rank % dev_count ;

  KokkosArray::Cuda::SelectDevice select_device( dev_rank );
  KokkosArray::Cuda::initialize( select_device );
  Explicit::driver<double,KokkosArray::Cuda>( "Cuda" , machine , elem_count_begin , elem_count_end , uq_count_begin , uq_count_end , count_run );
  KokkosArray::Cuda::finalize();
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

