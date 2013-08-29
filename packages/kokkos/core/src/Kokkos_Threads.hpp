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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_THREADS_HPP
#define KOKKOS_THREADS_HPP

#include <cstddef>
#include <iosfwd>
#include <Kokkos_Layout.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <Kokkos_HostSpace.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class ThreadsExec ;
} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Device for a pool of Pthreads or C11 threads on a CPU. */
class Threads {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  typedef Threads                  device_type ;
  typedef Kokkos::HostSpace        memory_space ;
  typedef memory_space::size_type  size_type ;
  typedef Kokkos::LayoutRight      array_layout ;
  typedef Kokkos::Threads          host_mirror_device_type ;

  //@}
  /*------------------------------------------------------------------------*/
  //! \name Static functions that all Kokkos devices must implement.
  //@{

  /** \brief  Query if called within a thread-parallel function */
  static int in_parallel();

  /** \brief  Set the device in a "sleep" state.
   *
   * This function sets the device in a "sleep" state in which it is
   * not ready for work.  This may consume less resources than if the
   * device were in an "awake" state, but it may also take time to
   * bring the device from a sleep state to be ready for work.
   *
   * \return True if the device is in the "sleep" state, else false if
   *   the device is actively working and could not enter the "sleep"
   *   state.
   */
  static bool sleep();

  /// \brief Wake the device from the 'sleep' state so it is ready for work.
  ///
  /// \return True if the device is in the "ready" state, else "false"
  ///  if the device is actively working (which also means that it's
  ///  awake).
  static bool wake();

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void fence();

  /// \brief Free any resources being consumed by the device.
  ///
  /// For the Threads device, this terminates spawned worker threads.
  static void finalize();

  /** \brief  Print configuration information */
  static void print_configuration( std::ostream & , bool detail = false );

  //@}
  /*------------------------------------------------------------------------*/
  /** \name Function for the functor device interface */
  /**@{ */

  inline int league_rank() const ;
  inline int league_size() const ;
  inline int team_rank() const ;
  inline int team_size() const ;

  inline void team_barrier();

  /* Collectively compute the league-wide unordered exclusive prefix sum.
   * Values are ordered within a team, but not between teams (i.e. the start
   * values of thread 0 in each team are not ordered according to team number).
   * This call does not use a global synchronization. Multiple unordered scans
   * can be in-flight at the same time (using different scratch_views).
   * The scratch-view will hold the complete sum in the end.
   */
  template< class VT >
  inline typename VT::value_type unordered_scan
             (typename VT::value_type& value, VT& scratch_view);

  /* Collectively compute the team-wide exclusive prefix sum using CUDA Unbound.
   * Values are ordered, the last thread returns the sum of all values
   * in the team less its own value
   */
  template< typename T >
  inline T team_scan(T& value);

  template< typename T >
  inline T * get_shmem( const int count );

  explicit inline Threads( Impl::ThreadsExec & );

  /**@} */
  /*------------------------------------------------------------------------*/
  //! \name Device-specific functions
  //@{

  /** \brief Initialize the device in the "ready to work" state.
   *
   *  The device is initialized in a "ready to work" or "awake" state.
   *  This state reduces latency and thus improves performance when
   *  dispatching work.  However, the "awake" state consumes resources
   *  even when no work is being done.  You may call sleep() to put
   *  the device in a "sleeping" state that does not consume as many
   *  resources, but it will take time (latency) to awaken the device
   *  again (via the wake()) method so that it is ready for work.
   *
   *  The 'team_topology' argument specifies
   *    (first) the number of teams of threads - the league_size
   *    (second) the number of threads per team - the team_size.
   *
   *  The 'use_core_topology' argument specifies the subset of available cores 
   *  to use for the device's threads.  If 'hwloc' is available then the
   *  full core topology can be queried via 'hwloc::get_core_topology()'.
   *  If 'use_core_topology' is not specified and 'hwloc' is available
   *  then the full core topology is used.  If 'hwloc' in not available
   *  then 'use_core_topology' is ignored.
   *
   *  Teams of threads are evenly distributed across the core topology.
   *  If team_topology.first*team_topology.second is less than or equal to
   *  use_core_topology.first*use_core_topology.second then each thread
   *  is assigned its own core.  Otherwise multiple threads in a team are
   *  assigned to a shared core (using hyperthreads).
   */
  static void initialize( const std::pair<unsigned,unsigned> team_topology ,
                          const std::pair<unsigned,unsigned> use_core_topology =
                                std::pair<unsigned,unsigned>(0u,0u) );


  static int league_max();
  static int team_max();

  //@}
  /*------------------------------------------------------------------------*/

private:

  friend class Impl::ThreadsExec ;

  Impl::ThreadsExec & m_exec ;
};

/*--------------------------------------------------------------------------*/

} // namespace Kokkos

#include <Kokkos_Parallel.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>
#include <Threads/Kokkos_Threads_Parallel.hpp>

#endif /* #define KOKKOS_THREADS_HPP */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

