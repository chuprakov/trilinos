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

#include <iostream>
#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void * host_allocate_not_thread_safe(
  const std::string    & label ,
  const std::type_info & scalar_type ,
  const size_t           scalar_size ,
  const size_t           scalar_count );

void host_decrement_not_thread_safe( const void * ptr );

}
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

enum { CURRENT_MEMORY_FUNCTOR_SIZE = 0x08000 /* 32k words , 128k bytes */ };

ThreadsExec                  s_threads_process ;
ThreadsExec                * s_threads_exec[  ThreadsExec::MAX_THREAD_COUNT ];
std::pair<unsigned,unsigned> s_threads_coord[ ThreadsExec::MAX_THREAD_COUNT ];
std::string                  s_exception_msg ;

unsigned s_threads_count       = 0 ;
unsigned s_threads_reduce_size = 0 ;
unsigned s_threads_shared_size = 0 ;

const void * s_current_function_lock = 0 ;

void (* volatile s_current_function)( ThreadsExec & , const void * );
const void * volatile s_current_function_arg = 0 ;

struct Sentinel {
  Sentinel() {}
  ~Sentinel()
  {
    if ( s_threads_count ||
         s_threads_reduce_size ||
         s_threads_shared_size ||
         s_current_function ||
         s_current_function_arg ||
         s_threads_exec[0] ) {
      std::cerr << "ERROR : Process exiting without calling Kokkos::Threads::terminate()" << std::endl ;
    }
  }
};

} // namespace
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsExec::driver(void)
{
  // If hardware locality library unavailable then pass in the rank.

  size_t thread_rank = (size_t) s_current_function_arg ;

  if ( s_threads_count <= thread_rank ) {
    thread_rank = Kokkos::hwloc::bind_this_thread( s_threads_count , s_threads_coord );
  }

  if ( s_threads_count <= thread_rank || 0 != s_threads_exec[ thread_rank ] ) {

    // An error occured. Inform process that thread is terminating
    s_threads_process.m_state = ThreadsExec::Terminating ;

    return ;
  }

  {
    ThreadsExec this_thread ;

    this_thread.m_state = ThreadsExec::Active ;

    s_threads_exec[ thread_rank ] = & this_thread ;

    // Inform spawning process that the threads_exec entry has been set.
    s_threads_process.m_state = ThreadsExec::Active ;

    while ( ThreadsExec::Active == this_thread.m_state ) {

#if 0
      try {
        // Call work function
        (*s_current_function)( this_thread , s_current_function_arg );
      }
      catch( const std::exception & x ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception : " << x.what() << std::endl ;
        s_exception_msg.append( msg.str() );
      }
      catch( ... ) {
        std::ostringstream msg ;
        msg << "Kokkos::Threads[" << thread_rank << "] Uncaught exeception"  << std::endl ;
        s_exception_msg.append( msg.str() );
      }
#else
        (*s_current_function)( this_thread , s_current_function_arg );
#endif

      // Deactivate thread and wait for reactivation
      this_thread.m_state = ThreadsExec::Inactive ;
      wait_yield( this_thread.m_state , ThreadsExec::Inactive );
    }

    s_threads_process.m_state = ThreadsExec::Terminating ;

    s_threads_exec[ thread_rank ] = 0 ;
  }
}

void execute_function_noop( ThreadsExec & , const void * ) {}

} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

ThreadsExec::~ThreadsExec()
{
  m_reduce      = 0 ;
  m_shared      = 0 ;
  m_shared_end  = 0 ;
  m_shared_iter = 0 ;
  m_state       = ThreadsExec::Terminating ;
  m_fan_size    = 0 ;
  m_fan_team_size = 0 ;
  m_team_rank   = 0 ;
  m_team_size   = 0 ;
  m_league_rank = 0 ;
  m_league_rank = 0 ;
  m_thread_size = 0 ;
  m_thread_size = 0 ;

  for ( unsigned i = 0 ; i < MAX_FAN_COUNT ; ++i ) { m_fan[i] = 0 ; }
}

ThreadsExec::ThreadsExec()
  : m_reduce(0)
  , m_shared(0)
  , m_shared_end(0)
  , m_shared_iter(0)
  , m_state( ThreadsExec::Terminating )
  , m_fan_size(0)
  , m_fan_team_size(0)
  , m_team_rank(0)
  , m_team_size(0)
  , m_league_rank(0)
  , m_league_size(0)
  , m_thread_rank(0)
  , m_thread_size(0)
{
  for ( unsigned i = 0 ; i < MAX_FAN_COUNT ; ++i ) { m_fan[i] = 0 ; }
}

void ThreadsExec::set_threads_relationships(
  const std::pair<unsigned,unsigned> team_topo ,
  ThreadsExec * threads[] )
{
  const unsigned thread_count = team_topo.first * team_topo.second ;

  for ( unsigned r = 0 ; r < thread_count ; ++r ) {
    if ( threads[r] == 0 ) {
      Kokkos::Impl::throw_runtime_exception( std::string("ThreadsExec::set_threads_relationships FAILED : NULL entry" ) );
    }
  }

  for ( unsigned league_rank = 0 , r = 0 ; league_rank < team_topo.first ;  ++league_rank ) {
  for ( unsigned team_rank = 0 ;           team_rank   < team_topo.second ; ++team_rank , ++r ) {

    ThreadsExec & th = * threads[r] ;

    th.m_league_rank = league_rank ;
    th.m_league_size = team_topo.first ;
    th.m_team_rank   = team_rank ;
    th.m_team_size   = team_topo.second ;
    th.m_thread_rank = th.m_team_rank + th.m_team_size * th.m_league_rank ;
    th.m_thread_size = th.m_team_size * th.m_league_size ;

    th.m_fan_size = 0 ;
    th.m_fan_team_size = 0 ;

    // Intra-team reduction:
    for ( int n = 1 ; ( th.m_team_rank + n < th.m_team_size ) &&
                      ( 0 == ( n & th.m_team_rank ) ) ; n <<= 1 ) {
      th.m_fan[ th.m_fan_size ] = threads[ ( th.m_team_rank + n ) + ( th.m_league_rank * th.m_team_size ) ];
      ++th.m_fan_size ;
      ++th.m_fan_team_size ;
    }

    // Inter-team (intra-league) reduction:

    if ( th.m_team_rank == 0 ) {

      for ( int n = 1 ; ( th.m_league_rank + n < th.m_league_size ) &&
                        ( 0 == ( n & th.m_league_rank ) ) ; n <<= 1 ) {

        th.m_fan[ th.m_fan_size++ ] = threads[ ( th.m_league_rank + n ) * th.m_team_size ];
      }
    }
  }}
}

void ThreadsExec::execute_get_binding( ThreadsExec & exec , const void * )
{
  s_threads_coord[ exec.m_thread_rank ] = Kokkos::hwloc::get_this_thread_coordinate();
}

void ThreadsExec::execute_sleep( ThreadsExec & exec , const void * )
{
  ThreadsExec::global_lock();
  ThreadsExec::global_unlock();

  const int n = exec.m_fan_size ;

  for ( int i = 0 ; i < n ; ++i ) {
    wait( exec.m_fan[i]->m_state , ThreadsExec::Active );
  }

  exec.m_state = ThreadsExec::Inactive ;
}

void ThreadsExec::execute_reduce_resize( ThreadsExec & exec , const void * )
{
  if ( exec.m_reduce ) {
    Kokkos::Impl::host_decrement_not_thread_safe( exec.m_reduce );
    exec.m_reduce = 0 ;
  }

  if ( s_threads_reduce_size ) {

    exec.m_reduce =
      Kokkos::Impl::host_allocate_not_thread_safe( "reduce_scratch_space" , typeid(unsigned char) , 1 , s_threads_reduce_size );

    // Guaranteed multiple of 'unsigned'

    unsigned * ptr = (unsigned *)( exec.m_reduce );
    unsigned * const end = ptr + s_threads_reduce_size / sizeof(unsigned);

    // touch on this thread
    while ( ptr < end ) *ptr++ = 0 ;
  }
}

void ThreadsExec::execute_shared_resize( ThreadsExec & exec , const void * )
{
  if ( exec.m_team_rank ) {
    exec.m_shared = 0 ;
  }
  else {

    if ( exec.m_shared ) {
      Kokkos::Impl::host_decrement_not_thread_safe( exec.m_shared );
      exec.m_shared = 0 ;
    }

    if ( s_threads_shared_size ) {

      exec.m_shared =
        Kokkos::Impl::host_allocate_not_thread_safe( "shared_scratch_space" , typeid(unsigned char) , 1 , s_threads_shared_size );

      // Guaranteed multiple of 'unsigned'

      unsigned * ptr = (unsigned *)( exec.m_shared );
      unsigned * const end = ptr + s_threads_shared_size / sizeof(unsigned);

      // touch on this thread
      while ( ptr < end ) *ptr++ = 0 ;
    }
  }

  exec.m_shared_end = s_threads_shared_size ;
}

void * ThreadsExec::get_shmem( const int size )
{
  enum { MASK  = ParallelWorkRequest::shared_align - 1 };

  // m_shared_iter is in bytes, convert to integer offsets
  const int offset = m_shared_iter >> power_of_two<sizeof(int)>::value ;

  // Round up size alignment
  m_shared_iter += ( size + int(MASK) ) & ~int(MASK);

  if ( m_shared_end < m_shared_iter ) {
    Kokkos::Impl::throw_runtime_exception( std::string("ThreadsExec::get_shmem FAILED : exceeded shared memory size" ) );
  }

  return ((int*)m_shared) + offset ;
}

}
}

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

void ThreadsExec::verify_is_process( const std::string & name , const bool initialized )
{
  if ( ! is_process() ) {
    std::string msg( name );
    msg.append( " FAILED : Called by a worker thread, can only be called by the master process." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }

  if ( initialized && 0 == s_threads_count ) {
    std::string msg( name );
    msg.append( " FAILED : Threads not initialized." );
    Kokkos::Impl::throw_runtime_exception( msg );
  }
}

// Wait for root thread to become inactive
void ThreadsExec::fence()
{
  if ( s_threads_count ) {
    wait( s_threads_exec[0]->m_state , ThreadsExec::Active );

    s_current_function     = 0 ;
    s_current_function_arg = 0 ;

    if ( s_exception_msg.size() ) {
      Kokkos::Impl::throw_runtime_exception( s_exception_msg );
    }
  }
}

inline
void ThreadsExec::activate_threads()
{
  for ( unsigned i = s_threads_count ; 0 < i ; ) {
    --i ;
    s_threads_exec[i]->m_state = ThreadsExec::Active ;
  }
}

/** \brief  Begin execution of the asynchronous functor */
void ThreadsExec::start( void (*func)( ThreadsExec & , const void * ) , const void * arg )
{
  verify_is_process("ThreadsExec::start" , true );

  if ( s_current_function_lock && 
       s_current_function_lock != arg ) {
    Kokkos::Impl::throw_runtime_exception( std::string( "ThreadsExec::start() FAILED : otherwise acquired" ) );
  }

  s_exception_msg.clear();

  s_current_function     = func ;
  s_current_function_arg = arg ;

  activate_threads();

  if ( s_threads_process.m_thread_size ) {
    (*func)( s_threads_process , arg );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }
}

void * ThreadsExec::execute( void (*func)( ThreadsExec & , const void * ) , const void * arg )
{
  verify_is_process("ThreadsExec::execute" , true );

  if ( s_current_function || s_current_function_arg ) {
    Kokkos::Impl::throw_runtime_exception( std::string( "ThreadsExec::execute() FAILED : already executing" ) );
  }

  s_current_function     = func ;
  s_current_function_arg = arg ;

  for ( unsigned i = s_threads_count ; 0 < i ; ) {
    --i ;
    s_threads_exec[i]->m_state = ThreadsExec::Active ;
  }

  if ( s_threads_process.m_thread_size ) {
    (*func)( s_threads_process , arg );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }
  else {
    wait( s_threads_exec[0]->m_state , ThreadsExec::Active );
  }

  s_current_function     = 0 ;
  s_current_function_arg = 0 ;

  return s_threads_exec[0]->m_reduce ;
}

//----------------------------------------------------------------------------

bool ThreadsExec::sleep()
{
  verify_is_process("ThreadsExec::sleep", true );

  if ( & execute_sleep == s_current_function ) return false ;

  fence();

  ThreadsExec::global_lock();

  s_exception_msg.clear();

  s_current_function = & execute_sleep ;

  activate_threads();
  
  return true ;
}

bool ThreadsExec::wake()
{
  verify_is_process("ThreadsExec::wake", true );

  if ( & execute_sleep != s_current_function ) return false ;

  ThreadsExec::global_unlock();

  if ( s_threads_process.m_thread_size ) {
    execute_sleep( s_threads_process , 0 );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }

  fence();

  return true ;
}

//----------------------------------------------------------------------------

void ThreadsExec::execute_serial( void (*func)( ThreadsExec & , const void * ) )
{
  s_exception_msg.clear();

  s_current_function = func ;

  const unsigned begin = s_threads_process.m_thread_size ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i ; ) {
    ThreadsExec & th = * s_threads_exec[ --i ];

    th.m_state = ThreadsExec::Active ;

    wait_yield( th.m_state , ThreadsExec::Active );
  }

  if ( s_threads_process.m_thread_size ) {
    s_threads_process.m_state = ThreadsExec::Active ;
    (*func)( s_threads_process , 0 );
    s_threads_process.m_state = ThreadsExec::Inactive ;
  }

  s_current_function = 0 ;
}

//----------------------------------------------------------------------------

void * ThreadsExec::root_reduce_scratch()
{
  return s_threads_count ? s_threads_exec[0]->m_reduce : (void*) 0 ;
}

void ThreadsExec::resize_reduce_scratch( size_t size )
{
  fence();

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_threads_reduce_size < size || 0 == size ) {

    s_threads_reduce_size = size ;

    execute_serial( & execute_reduce_resize );
  }
}

void ThreadsExec::resize_shared_scratch( size_t size )
{
  fence();

  const size_t rem = size % Kokkos::Impl::MEMORY_ALIGNMENT ;

  if ( rem ) size += Kokkos::Impl::MEMORY_ALIGNMENT - rem ;

  if ( s_threads_shared_size < size || 0 == size ) {
    s_threads_shared_size = size ;

    execute_serial( & execute_shared_resize );
  }

  for ( unsigned i = 0 ; i < s_threads_count ; ) {
    ThreadsExec & team_th = * s_threads_exec[i] ;

    for ( int j = 0 ; j < team_th.m_team_size ; ++j , ++i ) {
      s_threads_exec[i]->m_shared = team_th.m_shared ;
    }
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::print_configuration( std::ostream & s , const bool detail )
{
  verify_is_process("ThreadsExec::print_configuration",false);

  fence();

  const std::pair<unsigned,unsigned> core_topo = Kokkos::hwloc::get_core_topology();
  const unsigned core_size = Kokkos::hwloc::get_core_capacity();

  s << "Kokkos::Threads hwloc[" << core_topo.first << "x" << core_topo.second << "x" << core_size << "]" ;

  if ( s_threads_exec[0] ) {
    s << " team_league[" << s_threads_exec[0]->m_league_size << "x" << s_threads_exec[0]->m_team_size << "]" ;
    if ( 0 == s_threads_process.m_thread_size ) { s << " Asynchronous" ; }
    s << " ReduceScratch[" << s_threads_reduce_size << "]"
      << " SharedScratch[" << s_threads_shared_size << "]" ;
    s << std::endl ;

    if ( detail ) {

      execute_serial( & execute_get_binding );

      for ( unsigned i = 0 ; i < s_threads_count ; ++i ) {
        ThreadsExec * const th = s_threads_exec[i] ;
        s << "  Thread hwloc("
          << s_threads_coord[i].first << ","
          << s_threads_coord[i].second << ")" ;

        s_threads_coord[i].first  = ~0u ;
        s_threads_coord[i].second = ~0u ;

        if ( th ) {
          s << " rank(" << th->m_league_rank << "." << th->m_team_rank << ")" ;
          if ( th->m_fan_size ) {
            s << " Fan ranks" ;
            for ( int j = 0 ; j < th->m_fan_size ; ++j ) {
              s << " (" << th->m_fan[j]->m_league_rank << "." << th->m_fan[j]->m_team_rank << ")" ;
            }
          }
        }
        s << std::endl ;
      }
    }
  }
  else {
    s << " not initialized" << std::endl ;
  }
}

//----------------------------------------------------------------------------

void ThreadsExec::initialize( 
  const std::pair<unsigned,unsigned> team_topo ,
        std::pair<unsigned,unsigned> core_use )
{
  static const Sentinel sentinel ;

  verify_is_process("ThreadsExec::initialize",false);

  std::ostringstream msg ;

  msg << "Kokkos::Threads::initialize("
      << " team_topo(" << team_topo.first << "," << team_topo.second << ")"
      << ", core_use(" << core_use.first << "," << core_use.second << ")"
      << " )" ;

  if ( s_threads_count ) {
    msg << " FAILED : Already initialized" ;
    Kokkos::Impl::throw_runtime_exception( msg.str() );
  }

  //------------------------------------
  // Query hardware topology and capacity, if available.

  const unsigned                     thread_count = team_topo.first * team_topo.second ;
  const bool                         hwloc_avail  = Kokkos::hwloc::available();
  const std::pair<unsigned,unsigned> core_topo    = Kokkos::hwloc::get_core_topology();
  const unsigned                     core_cap     = Kokkos::hwloc::get_core_capacity();
  const unsigned                     capacity     = hwloc_avail ? core_topo.first * core_topo.second * core_cap : 0 ;
        std::pair<unsigned,unsigned> master_coord = Kokkos::hwloc::get_this_thread_coordinate();
        bool                         asynchronous = false ;

  //------------------------------------
  // Use HWLOC to determine coordinates for pinning threads.

  if ( hwloc_avail ) {

    if ( capacity         < thread_count ||
         core_topo.first  < core_use.first ||
         core_topo.second < core_use.second ) {
      msg << " FAILED : Requested more cores or threads than HWLOC reports are available "
          << " core_topology(" << core_topo.first << "," << core_topo.second << ")"
          << " thread_capacity(" << capacity << ")" ;
      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }

    if ( 0 == core_use.first || 0 == core_use.second ) {
      // User requested that we determine best use of cores.

      // Start by assuming use of all available cores
      core_use = core_topo ;

      if ( thread_count <= ( core_topo.first - 1 ) * core_topo.second ) {
        // Can spawn all requested threads on their own (NUMA) group of cores,
        // can execute asynchronously.
        --core_use.first ;
      }
      else if ( thread_count <= core_topo.first * ( core_topo.second - 1 ) ) {
        // Can spawn all requested threads on their own core and have excess core,
        // can execute asynchronously.
        --core_use.second ;
      }
      else if ( core_topo.first * core_topo.second < thread_count &&
                thread_count <= core_topo.first * ( core_topo.second - 1 ) * core_cap ) {
        // Will oversubscribe cores and can omit one core
        --core_use.second ;
      }
    }

    if ( core_use.first < core_topo.first ) {
      // Can omit a (NUMA) group of cores and execute work asynchronously
      // on the other groups.

      Kokkos::Impl::host_thread_mapping( team_topo , core_use , core_topo , s_threads_coord );

      // Don't use master thread's first core coordinate (NUMA region).
      // Originally mapped:
      //   begin = core_topo.first - core_use.first ;
      //   end   = core_topo.first ;
      // So can decrement.

      for ( unsigned i = 0 ; i < thread_count ; ++i ) {
        if ( s_threads_coord[i].first <= master_coord.first ) {
          --( s_threads_coord[i].first );
        }
      }

      asynchronous = true ;
    }
    else if ( core_use.second < core_topo.second ) {
      // Can omit a core from each group and execute work asynchronously

      Kokkos::Impl::host_thread_mapping( team_topo , core_use , core_topo , s_threads_coord );

      // Force master thread onto the highest rank unused core.
      master_coord.second = ( core_topo.second - core_use.second ) - 1 ;

      asynchronous = true ;
    }
    else {
      // Spawn threads with root thread on the master process' core

      Kokkos::Impl::host_thread_mapping( team_topo , core_use , core_topo , master_coord , s_threads_coord );

      s_threads_coord[0] = std::pair<unsigned,unsigned>( ~0u , ~0u );
    }
  }

  //------------------------------------
  // Spawn threads

  {
    const unsigned thread_spawn_begin  = asynchronous ? 0 : 1 ;
    unsigned       thread_spawn_failed = 0 ;

    s_threads_count    = thread_count ;
    s_current_function = & execute_function_noop ; // Initialization work function

    // If not fully utilizing the capacity then spawn threads for asynchronous execution.

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {

      s_threads_process.m_state = ThreadsExec::Inactive ;

      // If hwloc available then spawned thread will choose its own rank,
      // otherwise specify the rank.
      s_current_function_arg = (void*)( hwloc_avail ? ~0u : i );

      // Spawn thread executing the 'driver()' function.
      // Wait until spawned thread has attempted to initialize.
      // If spawning and initialization is successfull then
      // an entry in 's_threads_exec' will be assigned.
      if ( ThreadsExec::spawn() ) {
        wait_yield( s_threads_process.m_state , ThreadsExec::Inactive );
      }
    }

    // Wait for all spawned threads to deactivate before zeroing the function.

    for ( unsigned i = thread_spawn_begin ; i < thread_count ; ++i ) {
      ThreadsExec * const th = s_threads_exec[i] ;
      if ( th ) {
        wait_yield( th->m_state , ThreadsExec::Active );
      }
      else {
        ++thread_spawn_failed ;
      }
    }

    s_current_function     = 0 ;
    s_current_function_arg = 0 ;

    if ( thread_spawn_failed ) {

      s_threads_count = 0 ;

      msg << " FAILED " << thread_spawn_failed << " attempts to spawn threads" ;

      Kokkos::Impl::throw_runtime_exception( msg.str() );
    }

    Kokkos::hwloc::bind_this_thread( master_coord );

    if ( thread_spawn_begin ) { // Include the master thread
      new( & s_threads_process ) ThreadsExec();
      s_threads_process.m_state = ThreadsExec::Inactive ;
      s_threads_exec[0] = & s_threads_process ;
    }
  }

  //------------------------------------
  // Initialize team topology and fan-in/out relationships:
  ThreadsExec::set_threads_relationships( team_topo , s_threads_exec );

  // Initial allocations:
  ThreadsExec::resize_reduce_scratch( 4096 );
  ThreadsExec::resize_shared_scratch( 4096 );
}

//----------------------------------------------------------------------------

void ThreadsExec::finalize()
{
  verify_is_process("ThreadsExec::finalize",false);

  fence();

  resize_reduce_scratch(0);
  resize_shared_scratch(0);

  const unsigned begin = s_threads_process.m_thread_size ? 1 : 0 ;

  for ( unsigned i = s_threads_count ; begin < i-- ; ) {

    if ( s_threads_exec[i] ) {

      s_threads_exec[i]->m_state = ThreadsExec::Terminating ;

      wait_yield( s_threads_process.m_state , ThreadsExec::Inactive );

      s_threads_process.m_state = ThreadsExec::Inactive ;
    }
  }

  if ( s_threads_process.m_thread_size ) {
    ( & s_threads_process )->~ThreadsExec();
    s_threads_exec[0] = 0 ;
  }

  Kokkos::hwloc::unbind_this_thread();

  s_threads_count = 0 ;
}

//----------------------------------------------------------------------------

} /* namespace Impl */
} /* namespace Kokkos */


