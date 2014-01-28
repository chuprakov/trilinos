
#include <utility>
#include <iostream>

#include <KokkosCore_config.h>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_Threads.hpp>
#if defined( KOKKOS_HAVE_OPENMP)
#include <Kokkos_OpenMP.hpp>
#endif
#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

#include <feint_fwd.hpp>

int main()
{
#if defined( KOKKOS_HAVE_PTHREAD )
  {
    // Use 4 cores per NUMA region, unless fewer available

    const unsigned use_numa_count     = Kokkos::hwloc::get_available_numa_count();
    const unsigned use_cores_per_numa = std::min( 4u , Kokkos::hwloc::get_available_cores_per_numa() );

    Kokkos::Threads::initialize( use_numa_count * use_cores_per_numa );

    std::cout << "feint< Threads , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , false >();

    std::cout << "feint< Threads , Usingtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Threads , true  >();

    Kokkos::Threads::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_OPENMP )
  {
    // Use 4 cores per NUMA region, unless fewer available

    const unsigned use_numa_count     = Kokkos::hwloc::get_available_numa_count();
    const unsigned use_cores_per_numa = std::min( 4u , Kokkos::hwloc::get_available_cores_per_numa() );

    Kokkos::OpenMP::initialize( use_numa_count * use_cores_per_numa );

    std::cout << "feint< Threads , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::OpenMP , false >();

    std::cout << "feint< Threads , Usingtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::OpenMP , true  >();

    Kokkos::OpenMP::finalize();
  }
#endif

#if defined( KOKKOS_HAVE_CUDA )
  {
    // Initialize Host mirror device
    Kokkos::Cuda::host_mirror_device_type::initialize(1);
    const unsigned device_count = Kokkos::Cuda::detect_device_count();

    // Use the last device:
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(device_count-1) );

    std::cout << "feint< Cuda , NotUsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , false >();

    std::cout << "feint< Cuda , UsingAtomic >" << std::endl ;
    Kokkos::Example::feint< Kokkos::Cuda , true  >();

    Kokkos::Cuda::finalize();
    Kokkos::Cuda::host_mirror_device_type::finalize();

  }
#endif
}

