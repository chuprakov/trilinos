/*
//@HEADER
// ************************************************************************
// 
//   KokkosArray: Manycore Performance-Portable Multidimensional Arrays
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

/*--------------------------------------------------------------------------*/

#include <iostream>
#include <limits>

/* KokkosArray interfaces */

#include <KokkosArray_Host.hpp>
#include <Host/KokkosArray_Host_Internal.hpp>

/*--------------------------------------------------------------------------*/
/* Third Party Libraries */

/* Hardware locality library: http://www.open-mpi.org/projects/hwloc/ */
#include <hwloc.h>

#define  REQUIRED_HWLOC_API_VERSION  0x000010300

#if HWLOC_API_VERSION < REQUIRED_HWLOC_API_VERSION
#error "Requires  http://www.open-mpi.org/projects/hwloc/  Version 1.3 or greater"
#endif

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {
namespace {

void print_bitmap( const hwloc_bitmap_t bitmap )
{
  std::cout << "{" ;
  for ( int i = hwloc_bitmap_first( bitmap ) ;
        -1 != i ; i = hwloc_bitmap_next( bitmap , i ) ) {
    std::cout << " " << i ;
  }
  std::cout << " }" ;
}

int all_cores_available( const hwloc_topology_t host_topology ,
                         const hwloc_obj_t      node ,
                         const hwloc_bitmap_t   proc_cpuset )
{
  if ( hwloc_bitmap_weight( proc_cpuset ) ) {

    const int core_count =
      hwloc_get_nbobjs_inside_cpuset_by_type( host_topology ,
                                            node->allowed_cpuset ,
                                            HWLOC_OBJ_CORE );

    for ( int j = 0 ; j < core_count ; ++j ) {

      const hwloc_obj_t core =
        hwloc_get_obj_inside_cpuset_by_type( host_topology ,
                                           node->allowed_cpuset ,
                                           HWLOC_OBJ_CORE , j );

#if 0
      if ( 0 == hwloc_bitmap_isincluded( core->allowed_cpuset ,
                                         proc_cpuset ) ) {
        return 0 ;
      }
#else

      // If process' cpuset intersects core's cpuset
      // then process can access this core.
      // Must use intersection instead of inclusion because
      // MPI may bind the process to only one of the core's hyperthreads.
      //
      // Assumption: if the process can access any hyperthread of the core
      // then it has ownership of the entire core.
      // This assumes that it would be performance-detrimental
      // to spawn more than one MPI process per core and use nested threading.
      if ( 0 == hwloc_bitmap_intersects( proc_cpuset ,
                                         core->allowed_cpuset ) ) {
        return 0 ;
      }
#endif
    }
  }

  return 1 ;
}

} // namespace
} /* namespace Impl */
} /* namespace KokkosArray */

/*--------------------------------------------------------------------------*/

namespace KokkosArray {
namespace Impl {

class HostInternalHWLOC : public HostInternal {
private:

  enum { MAX_NODE_COUNT = 1024 };

  hwloc_topology_t m_host_topology ;
  hwloc_obj_type_t m_node_type ;          // hwloc type for a "node"
  unsigned         m_node_core_count ;    // Cores per node
  unsigned         m_node_core_pu_count ; // Processing units per core per node
  unsigned         m_node_rank[ MAX_NODE_COUNT ];

  void thread_mapping(
    const unsigned thread_rank ,
    unsigned & gang_rank ,
    unsigned & worker_rank ,
    unsigned & node_rank ,
    unsigned & core_rank ,
    unsigned & pu_rank ) const ;

public:

  ~HostInternalHWLOC();
  HostInternalHWLOC();

  bool bind_thread( const unsigned thread_rank ) const ;
  void print_configuration( std::ostream & ) const ;
};

//----------------------------------------------------------------------------

HostInternal & HostInternal::singleton()
{
  static HostInternalHWLOC self ; return self ;
}

//----------------------------------------------------------------------------

void HostInternalHWLOC::thread_mapping(
  const unsigned thread_rank ,
  unsigned & gang_rank ,
  unsigned & worker_rank ,
  unsigned & node_rank ,
  unsigned & core_rank ,
  unsigned & pu_rank ) const
{
  // Which node -> core -> processing unit

  gang_rank   = thread_rank / HostInternal::m_worker_count ;
  worker_rank = thread_rank % HostInternal::m_worker_count ;

  // How many cores will be used:

  const unsigned min_worker_per_core =
    HostInternal::m_worker_count / m_node_core_count ;

  const unsigned max_worker_per_core =
    min_worker_per_core +
    ( HostInternal::m_worker_count % m_node_core_count ? 1 : 0 );

  const unsigned core_base =
    m_node_core_count * max_worker_per_core - HostInternal::m_worker_count ;

  const unsigned core_base_worker_count = core_base * min_worker_per_core ;

  // Use the upper range of accessible nodes.
  // Assumes that the lower range is preferred by the operating system.

  const unsigned gang_unused =
    ( HostInternal::m_gang_capacity - HostInternal::m_gang_count ); 

  node_rank = m_node_rank[ gang_rank + gang_unused ];

  // Which core:

  core_rank = 
    worker_rank < core_base_worker_count ?
    worker_rank / min_worker_per_core :
    core_base +
    ( worker_rank - core_base_worker_count ) / max_worker_per_core ;

  pu_rank = worker_rank % m_node_core_pu_count ;
}

void HostInternalHWLOC::print_configuration( std::ostream & s ) const
{
  static const char name_numa[] = "NUMA" ;
  static const char name_socket[] = "SOCKET" ;
  static const char name_core[] = "CORE" ;
  static const char name_pu[] = "PU" ;
  static const char name_unknown[] = "unknown" ;

  enum { OBJ_TYPE_COUNT = 4 };

  static const char * const obj_names[] = { name_numa , name_socket , name_core , name_pu };

  static const char * const name_node = 
    HWLOC_OBJ_NODE   == m_node_type ? name_numa : (
    HWLOC_OBJ_SOCKET == m_node_type ? name_socket : (
    HWLOC_OBJ_CORE   == m_node_type ? name_core : (
    HWLOC_OBJ_PU     == m_node_type ? name_pu : name_unknown )));

  static const hwloc_obj_type_t obj_types[] =
    { HWLOC_OBJ_NODE   /* NUMA region     */
    , HWLOC_OBJ_SOCKET /* hardware socket */
    , HWLOC_OBJ_CORE   /* hardware core   */
    , HWLOC_OBJ_PU     /* processing unit */
    };

  const int obj_counts[] =
    { hwloc_get_nbobjs_by_type( m_host_topology , HWLOC_OBJ_NODE )
    , hwloc_get_nbobjs_by_type( m_host_topology , HWLOC_OBJ_SOCKET )
    , hwloc_get_nbobjs_by_type( m_host_topology , HWLOC_OBJ_CORE )
    , hwloc_get_nbobjs_by_type( m_host_topology , HWLOC_OBJ_PU )
    };

  s << std::endl << "KokkosArray::Host HWLOC { " ;
  for ( int k = 0 ; k < OBJ_TYPE_COUNT ; ++k ) {
    s << " " << obj_names[k] << "(" << obj_counts[k] << ")" ;
  }
  s << std::endl ;

  hwloc_bitmap_t proc_cpuset = hwloc_bitmap_alloc();
  hwloc_bitmap_t test_cpuset = hwloc_bitmap_alloc();

  hwloc_get_cpubind( m_host_topology , proc_cpuset , HWLOC_CPUBIND_PROCESS );

  for ( int k = 0 ; k < OBJ_TYPE_COUNT ; ++k ) {

    const hwloc_obj_type_t type  = obj_types[k] ;
    const int              count = obj_counts[k];
    const char * const     name  = obj_names[k] ;

    s << "  " << name << " {" ;

    for ( int i = 0 ; i < count ; ++i ) {

      const hwloc_obj_t obj =
        hwloc_get_obj_by_type( m_host_topology , type , i );

      hwloc_bitmap_and( test_cpuset , proc_cpuset , obj->allowed_cpuset );

      const int obj_weight       = hwloc_bitmap_weight( obj->allowed_cpuset );
      const int intersect_weight = hwloc_bitmap_weight( test_cpuset );

      if ( intersect_weight ) {
        if ( intersect_weight == obj_weight ) {
          s << " " << i ;
        }
        else {
          s << " [ " << i << " : "
            << intersect_weight << "/" << obj_weight
            << " ]" ;
        }
      }
    }

    s << " }" << std::endl ;
  }
  s << "}" << std::endl ;

  hwloc_bitmap_free( proc_cpuset );
  hwloc_bitmap_free( test_cpuset );

  //------------------------------------

  if ( 1 < HostInternal::m_thread_count ) {
    s << "  thread binding {"
      << std::endl ;

    for ( Host::size_type i = 0 ; i < HostInternal::m_thread_count ; ++i ) {
      const HostThread & thread = * HostThread::get_thread(i);

      unsigned gang_rank = 0 ;
      unsigned worker_rank = 0 ;
      unsigned node_rank = 0 ;
      unsigned core_rank = 0 ;
      unsigned pu_rank = 0 ;

      thread_mapping( thread.rank() ,
                      gang_rank ,
                      worker_rank ,
                      node_rank ,
                      core_rank ,
                      pu_rank );

      s << "    { thread(" << thread.rank()
        << ") : gang(" << thread.gang_rank()
        << ") worker(" << thread.worker_rank()
        << ") } -> { " << name_node << "(" << node_rank
        << ")" ;

      if ( 1 < m_node_core_count ) {
        s << " CORE(" << core_rank << ")" ;
      }
      if ( 1 < m_node_core_pu_count ) {
        s << " PU(" << pu_rank << ")" ;
      }
      s << " }" << std::endl ;
    }
    s << "  }" << std::endl ;
  }
}

//----------------------------------------------------------------------------

bool HostInternalHWLOC::bind_thread( const unsigned thread_rank ) const
{
  bool result = true ;

  // Can only safely bind threads if topology was detected

  if ( HostInternal::m_gang_capacity ) {

    unsigned gang_rank = 0 ;
    unsigned worker_rank = 0 ;
    unsigned node_rank = 0 ;
    unsigned core_rank = 0 ;
    unsigned pu_rank = 0 ;

    thread_mapping( thread_rank ,
                    gang_rank ,
                    worker_rank ,
                    node_rank ,
                    core_rank ,
                    pu_rank );

    const hwloc_obj_t node =
      hwloc_get_obj_by_type( m_host_topology, m_node_type, node_rank );

    const hwloc_obj_t core =
      hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                           node->allowed_cpuset ,
                                           HWLOC_OBJ_CORE ,
                                           core_rank );

    const hwloc_obj_t pu =
      hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                           core->allowed_cpuset ,
                                           HWLOC_OBJ_PU ,
                                           pu_rank );

    result = 0 == hwloc_set_cpubind( m_host_topology ,
                                     pu->allowed_cpuset ,
                                     HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT );

    if ( result ) {

      hwloc_cpuset_t thread_cpuset = hwloc_bitmap_alloc();

      hwloc_get_cpubind( m_host_topology, thread_cpuset, HWLOC_CPUBIND_THREAD );

      result = hwloc_bitmap_isequal( thread_cpuset , pu->allowed_cpuset );

      hwloc_bitmap_free( thread_cpuset );
    }
  }

  return result ;
}

//----------------------------------------------------------------------------

HostInternalHWLOC::HostInternalHWLOC()
  : HostInternal()
  , m_host_topology()
  , m_node_type( HWLOC_OBJ_CORE )
  , m_node_core_count( 0 )
  , m_node_core_pu_count( 0 )
{
  hwloc_topology_init( & m_host_topology );
  hwloc_topology_load( m_host_topology );

  int node_count = 0 ;

  //------------------------------------------------------------------------
  // 'Node' level of hierarchical topology:
  {
    // Choose a hwloc type for a 'node' from, in search order, the following:
    static const hwloc_obj_type_t candidate_node_type[] =
      { HWLOC_OBJ_NODE   /* NUMA region     */
      , HWLOC_OBJ_SOCKET /* hardware socket */
      , HWLOC_OBJ_CORE   /* hardware core   */
      };

    enum { CANDIDATE_NODE_TYPE_COUNT =
             sizeof(candidate_node_type) / sizeof(hwloc_obj_type_t) };

    hwloc_bitmap_t proc_cpuset = hwloc_bitmap_alloc();

    hwloc_get_cpubind( m_host_topology , proc_cpuset , HWLOC_CPUBIND_PROCESS );

    for ( int k = 0 ;
          k < CANDIDATE_NODE_TYPE_COUNT && 0 == node_count ; ++k ) {

      const hwloc_obj_type_t type = candidate_node_type[k] ;

      const int count = hwloc_get_nbobjs_by_type( m_host_topology , type );

      if ( 1 < count ) {

        for ( int i = 0 ; i < count ; ++i ) {

          const hwloc_obj_t node =
            hwloc_get_obj_by_type( m_host_topology , type , i );

          if ( all_cores_available( m_host_topology , node , proc_cpuset ) ) {
            m_node_rank[ node_count++ ] = i ;
          }
        }
      }

      if ( node_count ) {
        m_node_type = type ;
      }
    }

    hwloc_bitmap_free( proc_cpuset );
  }

  if ( 0 == node_count ) return ; // Failed to detect 'node' type

  //------------------------------------------------------------------------
  // Subsequent levels of hierarchy:

  int node_symmetry = 1 ;
  int cache_line_size = -1 ;

  int core_per_node = -1 ;
  int pu_per_core = -1 ;

  for ( int i = 0 ; i < node_count ; ++i ) {

    const hwloc_obj_t node =
      hwloc_get_obj_by_type( m_host_topology , m_node_type , m_node_rank[i] );

    const int core_count =
      hwloc_get_nbobjs_inside_cpuset_by_type( m_host_topology ,
                                              node->allowed_cpuset ,
                                              HWLOC_OBJ_CORE );

    if ( -1 == core_per_node ) { core_per_node = core_count ; }

    if ( core_count != core_per_node ) { node_symmetry = false ; }

    if ( core_count < core_per_node ) core_per_node = core_count ;

    for ( int j = 0 ; j < core_count ; ++j ) {

      const hwloc_obj_t core =
        hwloc_get_obj_inside_cpuset_by_type( m_host_topology ,
                                             node->allowed_cpuset ,
                                             HWLOC_OBJ_CORE , j );

      // Processing units (a.k.a., hyperthreads)

      const int pu_count =
        hwloc_get_nbobjs_inside_cpuset_by_type( m_host_topology ,
                                                core->allowed_cpuset ,
                                                HWLOC_OBJ_PU );

      if ( -1 == pu_per_core ) { pu_per_core = pu_count ; }

      if ( pu_count != pu_per_core ) { node_symmetry = 0 ; }

      if ( pu_count < pu_per_core ) { pu_per_core = pu_count ; }

      // Use the largest cache line size
      // assuming the largest will be a multiple of the smallest...

      const hwloc_obj_t core_cache_info =
        hwloc_get_shared_cache_covering_obj( m_host_topology , core );

      if ( core_cache_info && core_cache_info->attr ) {

        if ( -1 == cache_line_size ) {
          cache_line_size = core_cache_info->attr->cache.linesize ;
        }

        if ( cache_line_size != (int) core_cache_info->attr->cache.linesize ) {
          node_symmetry = 0 ;
        }

        if ( cache_line_size < (int) core_cache_info->attr->cache.linesize ) {
          cache_line_size = (int) core_cache_info->attr->cache.linesize ;
        }
      }
    }
  }

  m_node_core_count    = core_per_node ;
  m_node_core_pu_count = pu_per_core ;

  HostInternal::m_worker_capacity = core_per_node * pu_per_core ;
  HostInternal::m_gang_capacity   = node_count ;

  if ( 0 < cache_line_size ) {
    HostInternal::m_cache_line_size = cache_line_size ;
  }
}

HostInternalHWLOC::~HostInternalHWLOC()
{
  hwloc_topology_destroy( m_host_topology );
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

} /* namespace Impl */
} /* namespace KokkosArray */

