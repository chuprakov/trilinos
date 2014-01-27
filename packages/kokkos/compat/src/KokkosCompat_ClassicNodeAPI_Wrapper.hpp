#ifndef KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP
#define KOKKOSCOMPAT_CLASSICNODEAPI_WRAPPER_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Kokkos_View.hpp>
#ifdef KOKKOS_HAVE_CUDA
  #ifndef KERNEL_PREFIX
    #ifdef __CUDACC__
    #define KERNEL_PREFIX __host__ __device__
    #endif
  #endif
  #include <Kokkos_Cuda.hpp>
#endif
#ifdef KOKKOS_HAVE_OPENMP
#include <Kokkos_OpenMP.hpp>
#endif
#ifdef KOKKOS_HAVE_PTHREAD
#include <Kokkos_Threads.hpp>
#endif

#include <KokkosCompat_View.hpp>
/*namespace KokkosClassic {
      enum ReadWriteOption {
        ReadWrite, // < Indicates that the view may be safely read and written.
        WriteOnly  // < Indicates that the contents of the view are undefined until set on the host.
      };
}*/

namespace Kokkos {
  namespace Compat {
  /** \brief A default implementation of the Node memory architecture for a single memory space allocated by standard library calls.
      \ingroup kokkos_node_api
   */
  template<class DeviceType>
  class KokkosDeviceWrapperNode {
    public:

      //! Indicates that parallel buffers allocated by this node are available for use on the host thread.
      typedef DeviceType device_type;

      static const bool isHostNode = true;
      static const bool isCUDANode = false;

      static int count;

      KokkosDeviceWrapperNode(Teuchos::ParameterList &pl) {
        Teuchos::ParameterList params = getDefaultParameters();
        params.setParameters(pl);
        const int curNumThreads = params.get<int>("Num Threads");
        const int curNumTeams = params.get<int>("Num Teams");
        const int curDevice = params.get<int>("Device");
        int verboseInt = params.get<int>("Verbose");
        bool verbose = (verboseInt != 0);
        if (verbose) {
          std::cout << "DeviceWrapperNode initializing with \"numthreads\" = "
                    << curNumThreads << ", \"numteams\" = " << curNumTeams
                    << " \"device\" = " << curDevice << std::endl;
        }
        if(count==0)
          init (curNumTeams,curNumThreads,curDevice);
        count++;
      }

      KokkosDeviceWrapperNode() {
        Teuchos::ParameterList params = getDefaultParameters();
        const int curNumThreads = params.get<int>("Num Threads");
        const int curNumTeams = params.get<int>("Num Teams");
        const int curDevice = params.get<int>("Device");
        int verboseInt = params.get<int>("Verbose");
        bool verbose = (verboseInt != 0);
        if (verbose) {
          std::cout << "DeviceWrapperNode initializing with \"numthreads\" = "
              << curNumThreads << ", \"numteams\" = " << curNumTeams
              << " \"device\" = " << curDevice << std::endl;
        }
        if(count==0)
          init (curNumTeams,curNumThreads,curDevice);
        count++;
      };

      ~KokkosDeviceWrapperNode();

      Teuchos::ParameterList getDefaultParameters() {
        Teuchos::ParameterList params;
        params.set("Verbose",     0);
        params.set("Num Threads", 1);
        params.set("Num Teams", 1);
        params.set("Device", 0);
        return params;
      }

      void init(int numteams, int numthreads, int device);

      template <class WDP>
      struct FunctorParallelFor {
        typedef DeviceType device_type;

        const WDP _c;
        const int _beg;

        FunctorParallelFor (const int beg, const WDP wd) :
          _c (wd), _beg (beg)
        {}

        KOKKOS_INLINE_FUNCTION
        void operator() (const int & i) const {
          _c.execute (i + _beg);
        }
      };

      template <class WDP>
      static void parallel_for(int beg, int end, WDP wd) {
        const FunctorParallelFor<WDP> f(beg,wd);
        int n = end-beg;
        Kokkos::parallel_for(n,f);
      }

      template <class WDP>
      struct FunctorParallelReduce {
        typedef DeviceType device_type;
        typedef typename WDP::ReductionType value_type;

        WDP _c;
        const int _beg;
        FunctorParallelReduce (const int beg, WDP wd) :
          _c (wd), _beg (beg) {}

        KOKKOS_INLINE_FUNCTION
        void operator() (const int & i, volatile value_type& value) const {
          value = _c.reduce(value, _c.generate(i+_beg));
        }

        KOKKOS_INLINE_FUNCTION
        void init( volatile value_type &update) const
        {
          update = _c.identity();
        }

        KOKKOS_INLINE_FUNCTION
        void join( volatile value_type &update ,
                          const volatile value_type &source ) const
        {
          update = _c.reduce(update, source);
        }

      };

      template <class WDP>
      static typename WDP::ReductionType
      parallel_reduce(int beg, int end, WDP wd) {
        typedef typename WDP::ReductionType ReductionType;
        ReductionType globalResult = wd.identity();
        const FunctorParallelReduce<WDP> f(beg,wd);
        int n = end-beg;
        Kokkos::parallel_reduce(n,f,globalResult);
        return globalResult;
      }

      inline void sync() const {DeviceType::fence();};
      //@{ Memory management

      /*! \brief Allocate a parallel buffer, returning it as a pointer ecnapsulated in an ArrayRCP.

          Dereferencing the returned ArrayRCP or its underlying pointer in general results in undefined
          behavior outside of parallel computations.

          The buffer will be automatically freed by the Node when no more references remain.

          @tparam T The data type of the allocate buffer. This is used to perform alignment and determine the number of bytes to allocate.
          @param[in] size The size requested for the parallel buffer, greater than zero.

          \post The method will return an ArrayRCP encapsulating a pointer. The underlying pointer may be used in parallel computation routines,
                and is guaranteed to have size large enough to reference \c size number of entries of type \c T.
      */
      template <class T> inline
      Teuchos::ArrayRCP<T> allocBuffer(size_t size) {
        Teuchos::ArrayRCP<T> buff;
        if (size > 0) {
          Kokkos::View<T*,DeviceType> view("ANodeBuffer",size);
          buff = Teuchos::arcp(view.ptr_on_device(), 0, size,
                               deallocator(view), false);
        }
        if (isHostNode == false) {
          MARK_COMPUTE_BUFFER(buff);
        }
        return buff;
      }

      /*! \brief Copy data to host memory from a parallel buffer.

          @param[in] size       The number of entries to copy from \c buffSrc to \c hostDest.
          @param[in] buffSrc    The parallel buffer from which to copy.
          @param[out] hostDest  The location in host memory where the data from \c buffSrc is copied to.

          \pre  \c size is non-negative.
          \pre  \c buffSrc has length at least <tt>size</tt>.
          \pre  \c hostDest has length equal to \c size.
          \post On return, entries in the range <tt>[0 , size)</tt> of \c buffSrc have been copied to \c hostDest entries in the range <tt>[0 , size)</tt>.
      */
      template <class T> inline
      void copyFromBuffer(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayView<T> &hostDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffSrc);
        }
        Teuchos::ArrayRCP<T> buffDest = Teuchos::arcpFromArrayView(hostDest);
        copyBuffers(size,buffSrc,buffDest);
      }

      /*! \brief Copy data to host memory from a parallel buffer.

          @param[in]  size        The number of entries to copy from \c hostSrc to \c buffDest.
          @param[in]  hostSrc     The location in host memory from where the data is copied.
          @param[out] buffDest    The parallel buffer to which the data is copied.

          \pre  \c size is non-negative.
          \pre  \c hostSrc has length equal to \c size.
          \pre  \c buffSrc has length at least <tt>size</tt>.
          \post On return, entries in the range <tt>[0 , size)</tt> of \c hostSrc are allowed to be written to. The data is guaranteed to be present in \c buffDest before it is used in a parallel computation.
      */
      template <class T> inline
      void copyToBuffer(size_t size, const Teuchos::ArrayView<const T> &hostSrc, const Teuchos::ArrayRCP<T> &buffDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffDest);
        }
        Teuchos::ArrayRCP<const T> buffSrc = Teuchos::arcpFromArrayView(hostSrc);
        copyBuffers<T>(size,buffSrc,buffDest);
      }

      /*! \brief Copy data between buffers.

        @param[in]     size     The size of the copy, greater than zero.
        @param[in]     buffSrc  The source buffer, with length at least as large as \c size.
        @param[in,out] buffDest The destination buffer, with length at least as large as \c size.

        \post The data is guaranteed to have been copied before any other usage of buffSrc or buffDest occurs.
      */
      template <class T> inline
      void copyBuffers(size_t size, const Teuchos::ArrayRCP<const T> &buffSrc, const Teuchos::ArrayRCP<T> &buffDest) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buffSrc);
          CHECK_COMPUTE_BUFFER(buffDest);
        }
        Teuchos::ArrayView<const T> av_src = buffSrc(0,size);
        Teuchos::ArrayView<T>       av_dst = buffDest(0,size);
        std::copy(av_src.begin(),av_src.end(),av_dst.begin());
      }

      //! \brief Return a const view of a buffer for use on the host.
      template <class T> inline
      Teuchos::ArrayRCP<const T> viewBuffer(size_t size, Teuchos::ArrayRCP<const T> buff) {
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buff);
        }
        return buff.persistingView(0,size);
      }

      /// \brief Return a non-const view of a buffer for use on the host.
      ///
      /// \param rw [in] 0 if read-only, 1 if read-write.  This is an int and not a KokkosClassic::ReadWriteOption, in order to avoid a circular dependency between KokkosCompat and KokkosClassic.
      template <class T> inline
      Teuchos::ArrayRCP<T> viewBufferNonConst (const int rw, size_t size, const Teuchos::ArrayRCP<T> &buff) {
        (void) rw; // Silence "unused parameter" compiler warning
        if (isHostNode == false) {
          CHECK_COMPUTE_BUFFER(buff);
        }
        return buff.persistingView(0,size);
      }

      inline void readyBuffers(Teuchos::ArrayView<Teuchos::ArrayRCP<const char> > buffers, Teuchos::ArrayView<Teuchos::ArrayRCP<char> > ncBuffers) {
#ifdef HAVE_KOKKOSCLASSIC_DEBUG
        if (isHostNode == false) {
          for (size_t i=0; i < (size_t)buffers.size(); ++i) {
            CHECK_COMPUTE_BUFFER(buffers[i]);
          }
          for (size_t i=0; i < (size_t)ncBuffers.size(); ++i) {
            CHECK_COMPUTE_BUFFER(ncBuffers[i]);
          }
        }
#endif
        (void)buffers;
        (void)ncBuffers;
      }


      //@}
  };

  #ifdef KOKKOS_HAVE_CUDA
    typedef  KokkosDeviceWrapperNode<Kokkos::Cuda> KokkosCudaWrapperNode;
  #endif

  #ifdef KOKKOS_HAVE_OPENMP
    typedef  KokkosDeviceWrapperNode<Kokkos::OpenMP> KokkosOpenMPWrapperNode;
  #endif

  #ifdef KOKKOS_HAVE_PTHREAD
    typedef  KokkosDeviceWrapperNode<Kokkos::Threads> KokkosThreadsWrapperNode;
  #endif
  }
} // end of namespace Kokkos
#endif
