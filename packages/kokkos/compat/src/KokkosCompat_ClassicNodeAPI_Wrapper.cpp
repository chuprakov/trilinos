#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>
#include <Kokkos_HostSpace.hpp>

namespace Kokkos {
  namespace Compat {

    // mfh 01 Jan 2014: These definitions of the class variable count
    // need to be inside the namespace.  Declaring them as "template<>
    // int Kokkos::Compat::KokkosCudaWrapperNode::count = 0" in the
    // global namespace is a C++11 extension and results in compiler
    // warnings with Clang 3.2 on MacOS X.
#ifdef KOKKOS_HAVE_CUDA
    template<> int KokkosCudaWrapperNode::count = 0;
#endif
#ifdef KOKKOS_HAVE_OPENMP
    template<> int KokkosOpenMPWrapperNode::count = 0;
#endif
#ifdef KOKKOS_HAVE_PTHREAD
    template<> int KokkosThreadsWrapperNode::count = 0;
#endif
    template<> int KokkosSerialWrapperNode::count = 0;

#ifdef KOKKOS_HAVE_PTHREAD
    template<>
    KokkosDeviceWrapperNode<Kokkos::Threads>::
    ~KokkosDeviceWrapperNode<Kokkos::Threads> ()
    {
      count--;
      if (count == 0 && Threads::is_initialized ()) {
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::Threads,Cuda::host_mirror_device_type>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif
          Threads::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Threads>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::Threads::is_initialized ()) {
        if(NumNUMA>0 && NumCoresPerNUMA>0)
          Kokkos::Threads::initialize ( NumThreads, NumNUMA, NumCoresPerNUMA );
        else if (NumNUMA > 0)
          Kokkos::Threads::initialize ( NumThreads, NumNUMA );
        else
          Kokkos::Threads::initialize ( NumThreads );
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Threads>::name () {
      return "Threads/Wrapper";
    }
#endif

#ifdef KOKKOS_HAVE_OPENMP
    template<>
    KokkosDeviceWrapperNode<Kokkos::OpenMP>::~KokkosDeviceWrapperNode<Kokkos::OpenMP>() {
      count--;
      if (count == 0 && OpenMP::is_initialized ()) {
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::OpenMP, Cuda::host_mirror_device_type>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif
        OpenMP::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::OpenMP>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::OpenMP::is_initialized ()) {
        if(NumNUMA>0 && NumCoresPerNUMA>0)
          Kokkos::OpenMP::initialize ( NumThreads, NumNUMA, NumCoresPerNUMA );
        else if (NumNUMA > 0)
          Kokkos::OpenMP::initialize ( NumThreads, NumNUMA );
        else
          Kokkos::OpenMP::initialize ( NumThreads );
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::OpenMP>::name () {
      return "OpenMP/Wrapper";
    }
#endif

    template<>
    KokkosDeviceWrapperNode<Kokkos::Serial>::~KokkosDeviceWrapperNode<Kokkos::Serial>() {
      count--;
      if (count == 0 && Serial::is_initialized ()) {
#ifdef KOKKOS_HAVE_CUDA
        if (! Impl::is_same<Kokkos::Serial, Cuda::host_mirror_device_type>::value ||
            KokkosDeviceWrapperNode<Kokkos::Cuda>::count == 0)
#endif
        //Don't try to kill me if HostSpace was already destroyed.
        //Typical reason: static global instance of node is used, which might get destroyed after
        //the static HostSpace is destroyed.
        if(Kokkos::NEVEREVERUSEMEIWILLFINDYOU::host_space_singleton_wrapper().size()>0)
          Serial::finalize ();
      }
    }

    template<>
    void KokkosDeviceWrapperNode<Kokkos::Serial>::
    init (int NumThreads, int NumNUMA, int NumCoresPerNUMA, int Device) {
      if (! Kokkos::Serial::is_initialized ()) {
          Kokkos::Serial::initialize ();
      }
    }

    template<>
    std::string KokkosDeviceWrapperNode<Kokkos::Serial>::name () {
      return "Serial/Wrapper";
    }


  }
}



