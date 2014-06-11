/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_BulkData_hpp
#define stk_mesh_BulkData_hpp

//----------------------------------------------------------------------

#include <stddef.h>                     // for size_t, NULL
#include <stdint.h>                     // for uint16_t
#include <algorithm>                    // for max
#include <functional>                   // for less, equal_to
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <list>                         // for list
#include <map>                          // for map, map<>::value_compare
#include <stk_mesh/base/Entity.hpp>     // for Entity, etc
#include <stk_mesh/base/EntityCommDatabase.hpp>  // for EntityCommDatabase
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Selector.hpp>   // for Selector
#include <stk_mesh/base/Trace.hpp>      // for TraceIfWatching, etc
#include <stk_mesh/base/Types.hpp>      // for MeshIndex, EntityRank, etc
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/EntityRepository.hpp>  // for EntityRepository, etc
#include <stk_util/parallel/DistributedIndex.hpp>  // for DistributedIndex
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/util/TrackingAllocator.hpp>  // for tracking_allocator, etc
#include <string>                       // for char_traits, string
#include <utility>                      // for pair
#include <vector>                       // for vector
#include "Shards_CellTopologyData.h"    // for CellTopologyData, etc
#include "boost/functional/hash/extensions.hpp"  // for hash
#include "boost/tuple/detail/tuple_basic.hpp"  // for get
#include "boost/unordered/unordered_map.hpp"  // for unordered_map
#include "mpi.h"                        // for ompi_communicator_t
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::size_type, etc
#include "stk_mesh/base/BucketConnectivity.hpp"  // for BucketConnectivity
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey, hash_value
#include "stk_mesh/base/FieldDataManager.hpp"
#include "stk_mesh/base/Relation.hpp"   // for Relation, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowAssert, etc
namespace sierra { namespace Fmwk { class MeshBulkData; } }
namespace sierra { namespace Fmwk { class MeshObjSharedAttr; } }
namespace stk { namespace mesh { class FieldBase; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { namespace impl { class Partition; } } }
namespace stk { namespace mesh { struct ConnectivityMap; } }


//----------------------------------------------------------------------

// Use macro below to enable metric gathering for get_buckets memoization
//#define GATHER_GET_BUCKETS_METRICS

namespace stk {
namespace mesh {

namespace impl {

class EntityRepository;

}

parallel::DistributedIndex::KeySpanVector convert_entity_keys_to_spans( const MetaData & meta );

struct EntityCommListInfo
{
  EntityKey key;
  Entity    entity; // Might be invalid if entity has been deleted
  int  owner;
};

typedef TrackedVectorMetaFunc<EntityCommListInfo, EntityCommTag>::type EntityCommListInfoVector;

#ifdef __IBMCPP__
typedef std::vector<std::vector<FastMeshIndex> > VolatileFastSharedCommMapOneRank;
#else
typedef TrackedVectorMetaFunc<
  TrackedVectorMetaFunc<FastMeshIndex, VolatileFastSharedCommMapTag>::type,
  VolatileFastSharedCommMapTag>::type VolatileFastSharedCommMapOneRank;
#endif

inline
bool operator<(const EntityKey& key, const EntityCommListInfo& comm)
{ return key < comm.key; }

inline
bool operator<(const EntityCommListInfo& comm, const EntityKey& key)
{ return comm.key < key; }

inline
bool operator<(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return lhs.key < rhs.key; }

inline
bool operator==(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return lhs.key == rhs.key; }

inline
bool operator!=(const EntityCommListInfo& lhs, const EntityCommListInfo& rhs)
{ return !(lhs == rhs); }

struct IsInvalid
{
  bool operator()(const EntityCommListInfo& comm) const
  {
    return comm.key == EntityKey();
  }
};

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \brief  Manager for an integrated collection of
 *          \ref stk::mesh::Entity "entities",
 *          \ref stk::mesh::Relation "entity relations", and
 *          \ref stk::mesh::Bucket "buckets" of
 *          \ref stk_mesh_field_data "field data".
 *
 *  Bulk data should be distributed among all processors.
 */
class BulkData {

public:

  //Power users only.
  //Call this right after construction, before any field-data has been allocated.
  //If you call this method too late (after any field-data has been allocated, it will have no effect.
  //It turns off field-data updating so that movement of entities between buckets etc., as is done during
  //mesh-setup, will not cause corresponding churn of field-data.
  //Once the mesh is initialized with entities and relations, turn on field-data by calling the
  //method 'allocate_field_data'.
  void deactivate_field_updating();
  bool is_field_updating_active() const { return m_keep_fields_updated; }

  inline static BulkData & get( const Bucket & bucket);
  inline static BulkData & get( const Ghosting & ghost);
  inline static BulkData & get( const impl::BucketRepository & bucket_repo );

  typedef boost::unordered_map<EntityKey, size_t> GhostReuseMap;
  typedef std::map<std::pair<EntityRank, Selector>, std::pair<size_t, size_t> > SelectorCountMap;
#ifdef __IBMCPP__
  // The IBM compiler is easily confused by complex template types...
  typedef BucketVector                                                   TrackedBucketVector;
  typedef std::map<std::pair<EntityRank, Selector>, TrackedBucketVector> SelectorBucketMap;
  typedef std::vector<VolatileFastSharedCommMapOneRank>                  VolatileFastSharedCommMap;
#else
  typedef TrackedVectorMetaFunc<Bucket*, SelectorMapTag>::type  TrackedBucketVector;
  typedef std::map<std::pair<EntityRank, Selector>, TrackedBucketVector,
                   std::less<std::pair<EntityRank, Selector> >,
                   tracking_allocator<std::pair<std::pair<EntityRank, Selector>, TrackedBucketVector>, SelectorMapTag> > SelectorBucketMap;
  typedef TrackedVectorMetaFunc<VolatileFastSharedCommMapOneRank, VolatileFastSharedCommMapTag>::type VolatileFastSharedCommMap;

#endif

  enum BulkDataSyncState { MODIFIABLE = 1 , SYNCHRONIZED = 2 };

  ~BulkData();

  /** \brief  Construct mesh bulk data manager conformal to the given
   *          \ref stk::mesh::MetaData "meta data manager" and will
   *          distribute bulk data over the given parallel machine.
   *
   *  - The maximum number of entities per bucket may be supplied.
   *  - The bulk data is in the synchronized or "locked" state.
   */
  BulkData(   MetaData & mesh_meta_data
            , ParallelMachine parallel
#ifdef SIERRA_MIGRATION
            , bool add_fmwk_data = false
#endif
            , ConnectivityMap const* arg_connectivity_map = NULL
            , FieldDataManager *field_dataManager = NULL
            );

  //------------------------------------
  /** \brief  The meta data manager for this bulk data manager. */
  const MetaData & mesh_meta_data() const { return m_mesh_meta_data ; }
        MetaData & mesh_meta_data()       { return m_mesh_meta_data ; }

  /** \brief  The parallel machine */
  ParallelMachine parallel() const { return m_parallel_machine ; }

  /** \brief  Size of the parallel machine */
  int parallel_size()   const { return m_parallel_size ; }

  /** \brief  Rank of the parallel machine's local processor */
  int parallel_rank()   const { return m_parallel_rank ; }

  const ConnectivityMap & connectivity_map() const { return m_bucket_repository.connectivity_map(); }

  //------------------------------------
  /** \brief  Bulk data has two states:
   *          guaranteed to be parallel synchronized or
   *          modification in progress and may be parallel inconsistent.
   */
  BulkDataSyncState synchronized_state() const { return m_sync_state ; }

  enum modification_optimization {
      MOD_END_COMPRESS_AND_SORT
    , MOD_END_SORT
  };

  /** \brief  Count of the number of times that the bulk data has been
   *          parallel synchronized.  This count gets updated with
   *          each call to 'modification_end'.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  /** \brief  Begin a modification phase during which the mesh bulk data
   *          could become parallel inconsistent.  This is a parallel
   *          synchronous call.  The first time this method is called
   *          the mesh meta data is verified to be committed and
   *          parallel consistent.  An exception is thrown if this
   *          verification fails.
   *
   *  \return  True if transitioned out of the guaranteed
   *           parallel consistent state to the "ok to modify" state.
   *           False if already in this state.
   */
  bool modification_begin();


  /** \brief  Parallel synchronization of modifications and
   *          transition to the guaranteed parallel consistent state.
   *
   *  Parallel synchronization of accumulated local modifications
   *  is probably an expensive operation.  Operations include:
   *  - Determining ownership and sharing of created entities.
   *  - Synchronizing entity membership in parts for shared entities.
   *  - Refreshing the shared entities ghosting (e.g. aura).
   *  - Updating ghost entities that have change part membership.
   *  - Sorting buckets' entities for a well-defined ordering.
   *
   *  \return  True if transitioned from the "ok to modify" state.
   *           False if already already in this state.
   *
   *  \exception  If modification resolution errors occur then
   *              a parallel-consistent exception will be thrown.
   */
  bool modification_end( modification_optimization opt = MOD_END_SORT );


  /** If field-data was set to not stay in sync with buckets as the mesh was populated,
   * (by calling 'deactivate_field_updating' right after construction) this call
   * causes field-data to be allocated and field-data updating is re-activated.
   * If field-data was already allocated and staying in sync, then this call is a no-op.
   */
  void allocate_field_data();

  void verify_relations(const Bucket & bucket, Bucket::size_type bucket_ordinal, EntityRank rank) const;

  bool final_modification_end()
  {
    const bool mod_flag =  modification_end();

    //call modification_begin and end one last time to free deleted entities
    modification_begin();
    modification_end();

    m_mesh_finalized = true;

    return mod_flag;
  }

  /** \brief  Give away ownership of entities to other parallel processes.
   *
   *  A parallel-synchronous operation while the mesh is in the
   *  ok-to-modify state.
   *
   *  Each owning process inputs a list of entities and the
   *  new owning process.  Upon completion of the call the owning
   *  processes have the newly owned entities as well as the closure
   *  of those entities (without changing the ownership of the closure
   *  entities).  If a previous owner no longer needs a
   *  changed-owner entity to support the closure of a still-owned
   *  entity then the changed-owner entity is deleted from that process.
   *  All ghosts of all entities effected by the changed ownerships
   *  deleted.
   */
  void change_entity_owner( const std::vector<EntityProc> & arg_change);

  /** \brief  Rotate the field data of multistate fields.
   *
   *  <PRE>
   *  Rotation of states:
   *    StateN   <- StateNP1 (StateOld <- StateNew)
   *    StateNM1 <- StateN   (StateNM1 <- StateOld)
   *    StateNM2 <- StateNM1
   *    StateNM3 <- StateNM2
   *    StateNM3 <- StateNM2
   *  </PRE>
   */
  void update_field_data_states();

  /** \brief  Copy field data from src entity to Dest entity
   *           - Fields that exist on the src that don't exist on the dest will
   *             be ignored
   *           - Fields that exist on the dest that don't exist on the src will
   *             be zeroed or initialized with the Field-specified initial-value.
   */
  void copy_entity_fields( Entity src, Entity dst)
  {
    //TODO fix const correctness for src
    MeshIndex & src_mesh_idx = mesh_index(src);
    MeshIndex & dst_mesh_idx = mesh_index(dst);

    //// Pre-upgrade stk_mesh did not have this restriction, and it was easy enough to remove.
    //    ThrowAssert(src_mesh_idx.bucket->entity_rank() == dst_mesh_idx.bucket->entity_rank());

    copy_entity_fields_callback(dst_mesh_idx.bucket->entity_rank(),
                         dst_mesh_idx.bucket->bucket_id(),
                         dst_mesh_idx.bucket_ordinal,
                         src_mesh_idx.bucket->bucket_id(),
                         src_mesh_idx.bucket_ordinal);
  }

  //------------------------------------
  /** \brief  Query all buckets of a given entity rank
   *  Don't call inside BucketRepository member functions!
   */
  const BucketVector & buckets( EntityRank rank ) const
  { return m_bucket_repository.buckets(rank); }

  typedef impl::EntityRepository::const_iterator const_entity_iterator;

  //iterator that traverses entities of the specified rank, in order of ascending global identifier
  const_entity_iterator begin_entities(EntityRank ent_rank) const
  {
    return m_entity_repo.begin_rank(ent_rank);
  }

  //end-iterator for entities of the specified rank
  const_entity_iterator end_entities(EntityRank ent_rank) const
  {
    return m_entity_repo.end_rank(ent_rank);
  }

  /** \brief  Get entity with a given rank and id */
  Entity get_entity( EntityRank ent_rank , EntityId entity_id ) const {
    if (!is_good_rank_and_id(ent_rank, entity_id)) {
        return Entity();
    }
    return m_entity_repo.get_entity( EntityKey(ent_rank, entity_id));
  }


  /** \brief  Get entity with a given key */
  Entity get_entity( const EntityKey key ) const  {
    return m_entity_repo.get_entity(key);
  }

  //------------------------------------
  /** \brief  Create or retrieve a locally owned entity of a
   *          given rank and id.
   *
   *  A parallel-local operation.
   *
   *  The entity is created as locally owned and a member of the input
   *  mesh parts.  The entity a member of the meta data's locally owned
   *  mesh part and the entity's owner_rank() == parallel_rank().
   *
   *  If two or more processes create an entity of the same rank
   *  and identifier then the sharing and ownership of these entities
   *  will be resolved by the call to 'modification_end'.
   */
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , const PartVector& parts);

  Entity declare_entity( EntityRank ent_rank , EntityId ent_id , Part& part)
  {
    PartVector parts(1, &part);
    return declare_entity( ent_rank, ent_id, parts);
  }

  /** This overloading of declare_entity that doesn't take a part
   * creates the new entity in the 'universal' part.
   */
  Entity declare_entity( EntityRank ent_rank , EntityId ent_id);


  void change_entity_id( EntityId id, Entity entity);

  /** \brief  Change the parallel-locally-owned entity's
   *          part membership by adding and/or removing parts
   *
   *  A parallel-local operation.
   *
   *  If the locally owned entity is shared or ghosted then
   *  the change will be propagated to the sharing or ghosting
   *  processes by modification_end.
   */
  void change_entity_parts( Entity entity,
      const PartVector & add_parts ,
      const PartVector & remove_parts = PartVector() )
  {
    change_entity_parts(entity,
                        add_parts.begin(), add_parts.end(),
                        remove_parts.begin(), remove_parts.end());
  }

  //Optional parameter 'always_propagate_internal_changes' is always true except when this function
  //is being called from the sierra-framework. The fmwk redundantly does its own propagation of the
  //internal part changes (mostly induced-part stuff), so it's a performance optimization to avoid
  //the propagation that stk-mesh does.
  void change_entity_parts( Entity entity,
                            PartVector::const_iterator begin_add_parts, PartVector::const_iterator end_add_parts,
                            PartVector::const_iterator begin_remove_parts, PartVector::const_iterator end_remove_parts,
                            bool always_propagate_internal_changes=true );

  /** \brief  Request the destruction an entity on the local process.
   *
   * \paragraph destroy_requirements  Requirements
   *
   *  An entity cannot be the 'to' member of a relation.
   *  These relations must first be explicitly removed or the
   *  'from' entity be explicitly destroyed.
   *
   * \paragraph destroy_locally_owned  Destroy Locally Owned
   *
   *  Destruction of entities in the 'locally_owned_part' schedules
   *  all ghost copies of that entity for destruction during
   *  modification_end.  If the entity is shared with
   *  another process and that process does not also destroy the
   *  entity then ownership of the entity will be transfered to
   *  a sharing process during modification_end.
   *
   * \paragraph destroy_globally_shared  Destroy Locally Used
   *
   *  Entities in the 'globally_shared_part' are deleted
   *  on the local process and removed from the sharing lists on
   *  other processes during modication_end.
   *
   * \paragraph destroy_ghosted  Destroy Ghosted
   *
   *  Entities not in the 'locally_owned_part' and 'globally_shared_part'
   *  are ghosted.
   *  These entities are removed from all ghosting lists
   *  during 'modification_end'.
   *
   *  \return  True if the request for destruction is accepted; i.e.,
   *           if the entity is not the 'to' member of a relation.
   */
  bool destroy_entity( Entity entity, bool was_ghost = false );

  //------------------------------------

  /** \brief Generate a set of entites with globally unique id's
   *
   *  Each processor fills a request vector asking for a number of new
   *  entities of the given ranks.
   *
   *  ex. request = { 0, 4,  8}
   *  request 0 entites of rank 0, 4 entites of rank 1, and 8 entites
   *  of rank 2
   */
  void generate_new_entities(const std::vector<size_t>& requests,
      std::vector<Entity>& requested_entities);

  //------------------------------------
  /** \brief  Declare a relation and its converse between
   *          entities in the same mesh.
   *
   *  A parallel-local mesh modificaton operation.
   *
   *  This mapping ( e_from , local_id ) -> e_to  must be unique.
   *
   *  Relations between entities induces part membership as follows.
   *  1) If 'e_from' is a member of 'part' and
   *     part.primary_entity_rank() == entity_rank(e_from)
   *     then 'e_to' has induced membership in 'part'.
   *  2) If there exists a part relation 'part_rel' such that
   *     'e_from' is a member of part_rel.m_root and
   *     the entity relation conforms to the part relation
   *     then 'e_to' has induced membership in part_rel.m_target.
   *
   * Note that relation-declarations must be symmetric across all
   * sharers of the involved entities within a modification cycle.
   */
  void declare_relation( Entity e_from ,
      Entity e_to ,
      const RelationIdentifier local_id,
      Permutation permutation = static_cast<Permutation>(0));

  /** \brief  Declare a collection of relations by simply iterating
   *          the input and calling declare_relation on each entry.
   */
  void declare_relation( Entity entity, const std::vector<Relation> & rel);

  /** \brief  Remove all relations between two entities.
   *
   *  If the relation induced a part membership for 'e_to' and 'e_to'
   *  is not shared with another processor then that part membership
   *  is removed if and only if there does not exist another relation
   *  inducing the same part membership.
   *  If 'e_to' is shared then the check for removing the induced
   *  relatinship does not occur for that entity until the call to
   *  'modification_end'.
   *  The local_id arg is used to differentiate the case when there are
   *  multiple relationships between e_from and e_to.
   *
   *  Returns true if we were able to destroy the relation.
   */
  bool destroy_relation( Entity e_from ,
                         Entity e_to,
                         const RelationIdentifier local_id );


  // Check if entity has a specific relation to an entity of subcell_rank
  bool relation_exist( const Entity entity, EntityRank subcell_rank, RelationIdentifier subcell_id )
  {
    bool found = false;
    Entity const * rel_entity_it = bucket(entity).begin(bucket_ordinal(entity),subcell_rank);
    const unsigned num_rel = bucket(entity).num_connectivity(bucket_ordinal(entity),subcell_rank);
    ConnectivityOrdinal const * rel_ord_it = bucket(entity).begin_ordinals(bucket_ordinal(entity),subcell_rank);

    for (unsigned i=0 ; i < num_rel ; ++i) {
      if (rel_ord_it[i] == static_cast<ConnectivityOrdinal>(subcell_id) && is_valid(rel_entity_it[i])) {
        found = true;
        break;
      }
    }

    return found;
  }

  /** \brief  Determine the polarity of the local side,
   *          more efficient if the local_side_id is known.
   */
  bool element_side_polarity( const Entity elem ,
      const Entity side , unsigned local_side_id ) const
  {
    // 09/14/10:  TODO:  tscoffe:  Will this work in 1D??
    const bool is_side = entity_rank(side) != stk::topology::EDGE_RANK;
    const CellTopologyData * const elem_top = get_cell_topology( bucket(elem) ).getCellTopologyData();

    const unsigned side_count = ! elem_top ? 0 : (
        is_side ? elem_top->side_count
            : elem_top->edge_count );

    ThrowErrorMsgIf( elem_top == NULL,
        "For Element[" << identifier(elem) << "], element has no defined topology");

    ThrowErrorMsgIf( static_cast<unsigned>(side_count) <= local_side_id,
        "For Element[" << identifier(elem) << "], " <<
        "side: " << identifier(side) << ", " <<
        "local_side_id = " << local_side_id <<
        " ; unsupported local_side_id");

    const CellTopologyData * const side_top =
        is_side ? elem_top->side[ local_side_id ].topology
            : elem_top->edge[ local_side_id ].topology ;

    const unsigned * const side_map =
        is_side ? elem_top->side[ local_side_id ].node
            : elem_top->edge[ local_side_id ].node ;

    Entity const *elem_nodes = begin_nodes(elem);
    Entity const *side_nodes = begin_nodes(side);
    const unsigned n = side_top->node_count;
    bool good = false ;
    for ( unsigned i = 0 ; !good && i < n ; ++i ) {
        good = true;
        for ( unsigned j = 0; good && j < n ; ++j ) {
          good = side_nodes[(j+i)%n] == elem_nodes[ side_map[j] ];
        }
    }
    return good ;
  }

  //------------------------------------
  //------------------------------------
  /** \brief  All entities with communication information. */
  const EntityCommListInfoVector & comm_list() const
  { return m_entity_comm_list; }

  VolatileFastSharedCommMapOneRank const& volatile_fast_shared_comm_map(EntityRank rank) const
  {
    ThrowAssert(synchronized_state() == SYNCHRONIZED);
    ThrowAssertMsg(rank < stk::topology::ELEMENT_RANK, "Cannot shared entities of rank: " << rank);
    return m_volatile_fast_shared_comm_map[rank];
  }

  //------------------------------------
  /** \brief  Query the shared-entity aura.
   *          Is likely to be stale if ownership or sharing has changed
   *          and the 'modification_end' has not been called.
   */
  Ghosting & shared_aura() const { return * m_ghosting[1] ; }

  /** \brief Asymmetric parallel relations for owner-to-ghosted mesh entities.
   *
   *  - A collective parallel operation that must have the
   *    same name on all processors of this distributed mesh.
   */
  Ghosting & create_ghosting( const std::string & name );

  /** \brief  Change the members of a ghosting list on the sending processor.
   *
   *  - A collective parallel operation.
   *  - The ghosting must belong to this mesh.
   *  - Cannot change the 'shared_aura' in this manner.
   *  - Add locally owned entities to the input ghosting on the given
   *    destination processor.  The closure of the input entities
   *    will be ghosted.
   *  - Request removal of ghosted entities on the ghosting processor.
   *    This request will only be honored if the ghosted entity is
   *    not in the closure of another ghosted entity which will remain
   *    in or be added to this ghosting.
   */
  void change_ghosting( Ghosting & ghosts,
                        const std::vector<EntityProc> & add_send ,
                        const std::vector<EntityKey> & remove_receive = std::vector<EntityKey>());

  // Clear all ghosts for a particular ghosting.
  void destroy_ghosting( Ghosting& ghost_layer );

  /** \brief  Empty every single Ghosting.
   *          Same result, but more efficient than, calling
   *          change_ghosting to remove every single ghosted entity.
   */
  void destroy_all_ghosting();

  /** \brief  Query all ghostings */
  const std::vector<Ghosting*> & ghostings() const { return m_ghosting ; }

  /** \brief  Entity Comm functions that are now moved to BulkData
   */
  PairIterEntityComm entity_comm(const EntityKey & key) const { return m_entity_comm_map.comm(key); }
  PairIterEntityComm entity_comm_sharing(const EntityKey & key) const { return m_entity_comm_map.sharing(key); }
  PairIterEntityComm entity_comm(const EntityKey & key, const Ghosting & sub ) const { return m_entity_comm_map.comm(key,sub); }
  bool entity_comm_insert(Entity entity, const EntityCommInfo & val) { return m_entity_comm_map.insert(entity_key(entity), val, parallel_owner_rank(entity)); }
  bool entity_comm_erase(  const EntityKey & key, const EntityCommInfo & val) { return m_entity_comm_map.erase(key,val); }
  bool entity_comm_erase(  const EntityKey & key, const Ghosting & ghost) { return m_entity_comm_map.erase(key,ghost); }
  void entity_comm_clear_ghosting(const EntityKey & key ) { m_entity_comm_map.comm_clear_ghosting(key); }
  void entity_comm_clear(const EntityKey & key) { m_entity_comm_map.comm_clear(key); }
  int entity_comm_owner(const EntityKey & key) const;

  // Comm-related convenience methods

  bool in_shared(EntityKey key) const { return !entity_comm_sharing(key).empty(); }

  bool in_shared(EntityKey key, int proc) const;

  bool in_receive_ghost( EntityKey key ) const;

  bool in_receive_ghost( const Ghosting & ghost , EntityKey entity ) const;

  bool in_send_ghost( EntityKey key) const;

  bool in_send_ghost( EntityKey key , int proc ) const;

  bool in_ghost( const Ghosting & ghost , EntityKey key , int proc ) const;

  void comm_procs( EntityKey key, std::vector<int> & procs ) const; //shared and ghosted entities
  void comm_shared_procs( EntityKey key, std::vector<int> & procs ) const; // shared entities

  void shared_procs_intersection( std::vector<EntityKey> & keys, std::vector<int> & procs ) const;

  void comm_procs( const Ghosting & ghost ,
                   EntityKey key, std::vector<int> & procs ) const;

  //
  // Entity queries
  //

  bool in_index_range(Entity entity) const
  {
    return entity.local_offset() < m_entity_states.size();
  }

  bool is_valid(Entity entity) const
  {
    return (entity.local_offset() < m_entity_states.size()) && (m_entity_states[entity.local_offset()] != Deleted);
  }

  size_t count_relations(Entity entity) const;

  bool has_no_relations(Entity entity) const;

  void entity_setter_debug_check(Entity entity) const
  {
    // The 0-th local_offset is special, it represents the invalid, 0-initialized entity.
    // Client should never try to set properties on this entity even though it's in the index range.
    ThrowAssert(entity.local_offset() > 0);
  }

  void entity_getter_debug_check(Entity entity) const
  {
    ThrowAssertMsg(in_index_range(entity) , "Entity has out-of-bounds offset: " << entity.local_offset() << ", maximum offset is: " << m_entity_states.size() - 1);
  }

  //
  // Entity getters
  //

  inline const MeshIndex& mesh_index(Entity entity) const
  {
#ifndef NDEBUG
    entity_getter_debug_check(entity);
#endif

    return m_mesh_indexes[entity.local_offset()];
  }

  inline MeshIndex& mesh_index(Entity entity)
  {
#ifndef NDEBUG
    entity_setter_debug_check(entity); // setter check due to non-const
#endif

    return m_mesh_indexes[entity.local_offset()];
  }

  EntityId identifier(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return m_entity_keys[entity.local_offset()].id();
  }

  EntityRank entity_rank(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return m_entity_keys[entity.local_offset()].rank();
  }

  EntityKey entity_key(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return m_entity_keys[entity.local_offset()];
  }

  EntityState state(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return static_cast<EntityState>(m_entity_states[entity.local_offset()]);
  }

  size_t synchronized_count(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return m_entity_sync_counts[entity.local_offset()];
  }

  Bucket & bucket(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return *mesh_index(entity).bucket;
  }

  Bucket * bucket_ptr(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return mesh_index(entity).bucket;
  }

  Bucket::size_type bucket_ordinal(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return mesh_index(entity).bucket_ordinal;
  }

  int parallel_owner_rank(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return bucket(entity).parallel_owner_rank(bucket_ordinal(entity));
  }

  unsigned local_id(Entity entity) const
  {
    entity_getter_debug_check(entity);

    return m_local_ids[entity.local_offset()];
  }

#ifdef SIERRA_MIGRATION

  //this typedef for FmwkId must use the same type as the typedef
  //in the sierra-framework header framewk/mesh/Fmwk_Id.h
  //In other words, if you change this one, make the same change to the
  //one in Fmwk_Id.h.
  typedef int FmwkId; //must be a signed type -- fmwk uses negative values sometimes
  FmwkId global_id(Entity entity) const
  {
    ThrowAssertMsg(m_add_fmwk_data, "BulkData::global_id() only works under Framework mode");
    entity_getter_debug_check(entity);

    return m_fmwk_global_ids[entity.local_offset()];
  }

  const RelationVector& aux_relations(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    entity_setter_debug_check(entity); // setter check due to side effects

    if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
      m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
    }
    return *m_fmwk_aux_relations[entity.local_offset()];
  }

  RelationVector& aux_relations(Entity entity)
  {
    ThrowAssert(m_add_fmwk_data);
    entity_setter_debug_check(entity); // setter check due to side effects

    if (m_fmwk_aux_relations[entity.local_offset()] == NULL) {
      m_fmwk_aux_relations[entity.local_offset()] = new RelationVector();
    }
    return *m_fmwk_aux_relations[entity.local_offset()];
  }

  const sierra::Fmwk::MeshObjSharedAttr* get_shared_attr(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    entity_getter_debug_check(entity);

    return m_fmwk_shared_attrs[entity.local_offset()];
  }

  int get_connect_count(Entity entity) const
  {
    ThrowAssert(m_add_fmwk_data);
    entity_getter_debug_check(entity);

    return m_fmwk_connect_counts[entity.local_offset()];
  }

#endif

  //
  // Entity setters
  //

  void set_mesh_index(Entity entity, Bucket * in_bucket, Bucket::size_type ordinal )
  {
    // The trace statement forces this method to be defined after Entity
    TraceIfWatching("stk::mesh::BulkData::set_mesh_index", LOG_ENTITY, entity_key(entity));

    entity_setter_debug_check(entity);

    if (in_bucket != NULL) {
      ThrowAssertMsg(in_bucket->size() >= ordinal, "Detected bad bucket/ordinal.");
    }
    MeshIndex &mesh_idx = mesh_index(entity);
    mesh_idx.bucket = in_bucket;
    mesh_idx.bucket_ordinal = ordinal;
  }

  void set_entity_key(Entity entity, EntityKey key)
  {
    entity_setter_debug_check(entity);

    m_entity_keys[entity.local_offset()] = key;
  }

  void set_state(Entity entity, EntityState entity_state)
  {
    entity_setter_debug_check(entity);

    m_entity_states[entity.local_offset()] = static_cast<uint16_t>(entity_state);
  }

  bool set_parallel_owner_rank(Entity entity, int in_owner_rank)
  {
    TraceIfWatching("stk::mesh::BulkData::set_entity_owner_rank", LOG_ENTITY, entity_key(entity));
    DiagIfWatching(LOG_ENTITY, entity_key(entity), "new owner: " << in_owner_rank);

    entity_setter_debug_check(entity);

    int &rank = bucket(entity).m_owner_ranks[bucket_ordinal(entity)];
    if ( in_owner_rank != rank ) {
      rank = in_owner_rank;
      modified(entity);
      return true;
    }
    return false;
  }

  void set_synchronized_count(Entity entity, size_t sync_count)
  {
    entity_setter_debug_check(entity);

    m_entity_sync_counts[entity.local_offset()] = sync_count;
  }

  void set_local_id(Entity entity, unsigned id)
  {
    entity_setter_debug_check(entity);

    m_local_ids[entity.local_offset()] = id;
  }

#ifdef SIERRA_MIGRATION
  void set_global_id(Entity entity, int id)
  {
    ThrowAssert(m_add_fmwk_data);
    entity_setter_debug_check(entity);

    m_fmwk_global_ids[entity.local_offset()] = id;
  }

  template <typename SharedAttr>
  void set_shared_attr(Entity entity, SharedAttr* attr)
  {
    ThrowAssert(m_add_fmwk_data);
    entity_setter_debug_check(entity);

    m_fmwk_shared_attrs[entity.local_offset()] = attr;
  }

  void set_connect_count(Entity entity, int count)
  {
    ThrowAssert(m_add_fmwk_data);
    entity_setter_debug_check(entity);

    m_fmwk_connect_counts[entity.local_offset()] = count;
  }

  void set_relation_orientation(Entity from, Entity to, ConnectivityOrdinal to_ord, unsigned to_orientation);

  void reserve_relation(Entity entity, const unsigned num);
  void erase_and_clear_if_empty(Entity entity, RelationIterator rel_itr);
  void internal_verify_initialization_invariant(Entity entity);

#endif

  /**
   * Mark this entity as modified (only changes from Unchanged
   * to Modified). Propagates the modification to higher-ranking
   * entities related to this entity. In other words, based on our
   * modification model, all entities that have modified_entity in their
   * closure must also be marked as modified.
   */
  void modified(Entity entity);

  //
  // Connectivity getter methods. For each entity, you can get connected entities
  // of any rank (i.e. node, edge, face, element, constraint).
  // (Some of those connectivities are empty, for example an element may not
  // have any connected edges, depending whether edges exist in the mesh or not.)
  //
  // For each rank of connected entities, you can get the number of connected entities,
  // pointers to the connected entities themselves,
  // ordinals of those connected entities,
  // and permutations of the connected entities.
  //

  // Nodal connectivities are usually contiguous within a bucket. For example, an element-bucket
  // has all connected nodes for those elements allocated in contiguous memory. This should not be
  // assumed for other kinds of connectivity.

  //These connectivity getter methods are implemented by macros which are located further
  //down in this header. (Search for BEGIN_END_PAIR.)

  inline Entity const* begin(Entity entity, EntityRank rank) const;
  inline Entity const* begin_nodes(Entity entity) const;
  inline Entity const* begin_edges(Entity entity) const;
  inline Entity const* begin_faces(Entity entity) const;
  inline Entity const* begin_elements(Entity entity) const;

  // The ordinal of a connected entity is that entity's local index on the entity it
  // is connected to, as defined by standard exodus conventions. For example, the
  // connected nodes of a hex-8 element will have ordinals in the range 0 .. 7.

  // Connected entities are stored in order of ascending ordinal.
  // For cases like element-node connectivity, where an element always has all of its nodes,
  // the array of connected nodes can be indexed by ordinal.
  // For cases like element-face connectivity, where an element may not have all faces defined,
  // the array of connected faces can not always be indexed by ordinal even
  // though those faces are sorted by ordinal.
  // e.g., a hex-8 element may only have two of its six possible faces.

  inline ConnectivityOrdinal const* begin_ordinals(Entity entity, EntityRank rank) const;
  inline ConnectivityOrdinal const* begin_node_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_edge_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_face_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* begin_element_ordinals(Entity entity) const;

  // The permutation of a connected entity is an integer type which is used
  // to store the polarity and orientation of the entity.

  inline Permutation const* begin_permutations(Entity entity, EntityRank rank) const;
  inline Permutation const* begin_node_permutations(Entity entity) const;
  inline Permutation const* begin_edge_permutations(Entity entity) const;
  inline Permutation const* begin_face_permutations(Entity entity) const;
  inline Permutation const* begin_element_permutations(Entity entity) const;

  unsigned num_connectivity(Entity entity, EntityRank rank) const;

  inline unsigned num_nodes(Entity entity) const;
  inline unsigned num_edges(Entity entity) const;
  inline unsigned num_faces(Entity entity) const;
  inline unsigned num_elements(Entity entity) const;

  unsigned count_valid_connectivity(Entity entity, EntityRank rank) const;
  unsigned count_valid_connectivity(Entity entity) const;

  inline Entity const* end(Entity entity, EntityRank rank) const;
  inline Entity const* end_nodes(Entity entity) const;
  inline Entity const* end_edges(Entity entity) const;
  inline Entity const* end_faces(Entity entity) const;
  inline Entity const* end_elements(Entity entity) const;
  Entity const* end_constraints(Entity entity) const
  { return end(entity, stk::topology::CONSTRAINT_RANK); }
  inline ConnectivityOrdinal const* end_ordinals(Entity entity, EntityRank rank) const;
  inline ConnectivityOrdinal const* end_node_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_edge_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_face_ordinals(Entity entity) const;
  inline ConnectivityOrdinal const* end_element_ordinals(Entity entity) const;
  inline Permutation const* end_permutations(Entity entity, EntityRank rank) const;
  inline Permutation const* end_node_permutations(Entity entity) const;
  inline Permutation const* end_edge_permutations(Entity entity) const;
  inline Permutation const* end_face_permutations(Entity entity) const;
  inline Permutation const* end_element_permutations(Entity entity) const;


  // Return index (offset) of query ordinal if found, num_connectivity otherwise.
  unsigned find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const;

  bool has_permutation(Entity entity, EntityRank rank) const;

  bool owned_closure(Entity entity) const
  { return m_closure_count[entity.local_offset()] > static_cast<uint16_t>(0); }

  void get_selected_nodes(stk::mesh::Selector selector, stk::mesh::EntityVector& nodes);

  size_t total_field_data_footprint(const FieldBase &f, EntityRank rank) const
  {
    return m_bucket_repository.total_field_data_footprint(f, rank);
  }

  size_t total_field_data_footprint(EntityRank rank) const;

  //reserves space for a new entity, or reclaims space from a previously-deleted entity
  size_t generate_next_local_offset(size_t preferred_offset = 0);

#ifdef SIERRA_MIGRATION
  //strictly a transition aid!!! don't add new usage of this!
  void set_fmwk_bulk_data(const sierra::Fmwk::MeshBulkData* fmwk_bulk_ptr)
  {
    m_fmwk_bulk_ptr = fmwk_bulk_ptr;
  }

  //strictly a transition aid!!! don't add new usage of this!
  const sierra::Fmwk::MeshBulkData* get_fmwk_bulk_data() const
  {
    return m_fmwk_bulk_ptr;
  }

  RelationIterator internal_begin_relation(Entity entity, const Relation::RelationType relation_type) const
  {
    ThrowAssert(m_add_fmwk_data);
    if (impl::internal_is_handled_generically(relation_type)) {
      ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
      return RelationIterator();
    }
    else {
      return aux_relations(entity).begin();
    }
  }

  RelationIterator internal_end_relation(Entity entity, const Relation::RelationType relation_type) const
  {
    ThrowAssert(m_add_fmwk_data);
    if (impl::internal_is_handled_generically(relation_type)) {
      ThrowErrorMsg("stk::Mesh::BulkData::internal_begin_relation(..) requests native stk::mesh relation type");
      return RelationIterator();
    }
    else {
      return aux_relations(entity).end();
    }
  }

  void compress_relation_capacity(Entity entity)
  {
    RelationVector &rels = aux_relations(entity);
    RelationVector tmp(rels);
    tmp.swap(rels);
  }

  bool add_fmwk_data() const { return m_add_fmwk_data; }
#endif

  // Do not call!
  void internal_change_entity_key(EntityKey old_key, EntityKey new_key, Entity entity);

  // Do not call!  Just for a legacy test!
  impl::EntityRepository &get_entity_repository() { return m_entity_repo; }

  // Print all mesh info
  void dump_all_mesh_info(std::ostream& out = std::cout) const;

  // memoized version
  BucketVector const& get_buckets(EntityRank rank, Selector const& selector) const;

  // non-memoized version.
  void get_buckets(EntityRank rank, Selector const& selector, BucketVector & output_buckets) const;

  //
  //  Get entities of the specified rank that satisify the input selector.
  //  Note entities are returned in bucket order, though no particular order should be relied on
  //
  void get_entities(EntityRank rank, Selector const& selector, EntityVector& output_entities) const;


private:

  void update_deleted_entities_container();
  void addMeshEntities(const std::vector< stk::parallel::DistributedIndex::KeyTypeVector >& requested_key_types,
         const std::vector<Part*> &rem, const std::vector<Part*> &add, std::vector<Entity>& requested_entities);

#ifndef DOXYGEN_COMPILE

  // Forbidden
  BulkData();
  BulkData( const BulkData & );
  BulkData & operator = ( const BulkData & );

  //
  // Members
  //

  /** \brief  Parallel index for entity keys */
  parallel::DistributedIndex          m_entities_index;
  impl::EntityRepository              m_entity_repo;

  // Simply a list of data for entities that are being communicated
  EntityCommListInfoVector            m_entity_comm_list;

  // The full database of comm info for all communicated entities.
  EntityCommDatabase m_entity_comm_map;

  // Only works outside of modification cycles.
  // m_volatile_fast_shared_comm_map[entity_rank][parallel_rank] -> FastMeshIndexes
  //   Means that the entities represented by FastMeshIndexes are shared with proc parallel_rank
  VolatileFastSharedCommMap m_volatile_fast_shared_comm_map;

  std::vector<Ghosting*>              m_ghosting; /**< Aura is [1] */

  std::list<size_t, tracking_allocator<size_t, DeletedEntityTag> >     m_deleted_entities;
  std::list<size_t, tracking_allocator<size_t, DeletedEntityTag> >     m_deleted_entities_current_modification_cycle;

  GhostReuseMap m_ghost_reuse_map;

  // Other information:
  MetaData &         m_mesh_meta_data;
  ParallelMachine    m_parallel_machine;
  int                m_parallel_size;
  int                m_parallel_rank;
  size_t             m_sync_count;
  BulkDataSyncState  m_sync_state;
  bool               m_meta_data_verified;
  bool               m_mesh_finalized;
#ifdef SIERRA_MIGRATION
  bool                              m_add_fmwk_data; // flag that will add extra data to buckets to support fmwk
  const sierra::Fmwk::MeshBulkData* m_fmwk_bulk_ptr;

public:
  mutable bool       m_check_invalid_rels; // TODO REMOVE

private:
#endif
  int m_num_fields;
  bool m_keep_fields_updated;


  // Arrays of Entity member-data indexed by entity.local_offset():
  std::vector<MeshIndex> m_mesh_indexes;

  std::vector<EntityKey>   m_entity_keys;
  std::vector<uint16_t>    m_entity_states;
  std::vector<uint16_t>    m_closure_count;
  std::vector<size_t>      m_entity_sync_counts;
  std::vector<unsigned>    m_local_ids;

#ifdef SIERRA_MIGRATION
  // Extra data that fmwk needs to have on an entity. These vectors are indexed by local offset.

  mutable std::vector<RelationVector* > m_fmwk_aux_relations;   // Relations that can't be managed by STK such as PARENT/CHILD
  std::vector<FmwkId>                   m_fmwk_global_ids;
  std::vector<const sierra::Fmwk::MeshObjSharedAttr*>              m_fmwk_shared_attrs;
  std::vector<unsigned short>           m_fmwk_connect_counts;
#endif

 // There will be one of these per bucket

  // Outer index is (m_num_fields * entity rank) + field_ordinal, inner index
  // is bucket id, pair defines num bytes of data per entity and the
  // data for that field on that bucket


  // Outer index is rank, inner is bucket-id. This contains *all* field
  // data for a bucket.

//  ContiguousFieldDataManager m_default_field_data_manager;
  DefaultFieldDataManager m_default_field_data_manager;
  FieldDataManager *m_field_data_manager;

  // Memoize user bucket requests
  mutable SelectorBucketMap m_selector_to_buckets_map;
#ifdef GATHER_GET_BUCKETS_METRICS
  mutable SelectorCountMap m_selector_to_count_map;
  mutable size_t m_num_memoized_get_buckets_calls;
  mutable size_t m_num_non_memoized_get_buckets_calls;
  size_t m_num_buckets_inserted_in_cache;
  size_t m_num_buckets_removed_from_cache;
  size_t m_num_modifications;
#endif

  impl::BucketRepository              m_bucket_repository; // needs to be destructed first!

  //
  // Internal methods
  //

#ifdef GATHER_GET_BUCKETS_METRICS
  void gather_and_print_get_buckets_metrics() const;
#endif

  void gather_and_print_mesh_partitioning() const;

  // Field callbacks

  void new_bucket_callback(EntityRank rank, const PartVector& superset_parts, size_t capacity, Bucket* new_bucket);
  void new_bucket_caching(EntityRank rank, Bucket* new_bucket);

  //
  //  "fields" is an optional argument, if present copy only the listed fields.
  //
  void copy_entity_fields_callback(EntityRank dst_rank, unsigned dst_bucket_id, Bucket::size_type dst_bucket_ord,
                                   unsigned src_bucket_id, Bucket::size_type src_bucket_ord,
                                   const std::vector<FieldBase*>* fields =NULL);


  void destroy_bucket_callback(EntityRank rank, Bucket const& dying_bucket, unsigned capacity);

  // id_map, indexed by new id, maps to old id
  void reorder_buckets_callback(EntityRank rank, const std::vector<unsigned>& id_map);

  void remove_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);
  void remove_entity_field_data_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);
  void add_entity_callback(EntityRank rank, unsigned bucket_id, Bucket::size_type bucket_ord);

  // Misc

  void initialize_arrays();

  bool internal_declare_relation(Entity e_from, Entity e_to,
                                 RelationIdentifier local_id,
                                 unsigned sync_count, bool is_back_relation,
                                 Permutation permut);

  /** \brief  The meta data manager for this bulk data manager. */
  MetaData & meta_data() const { return m_mesh_meta_data ; }

  /** methods for managing arrays of entity member-data */

  void log_created_parallel_copy(Entity entity)
  {
    if (state(entity) == Created) {
      set_state(entity, Modified);
    }
  }

  /**
   * For all processors sharing an entity, find one to be the new
   * owner.
   */
  int determine_new_owner( Entity ) const ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */
  void internal_change_entity_parts( Entity ,
                                     const std::vector<Part*> & add_parts ,
                                     const std::vector<Part*> & remove_parts,
                                     bool always_propagate_internal_changes=true);

  void internal_propagate_part_changes( Entity entity, const std::vector<Part*> & removed );

  void internal_change_ghosting( Ghosting & ghosts,
                                 const std::vector<EntityProc> & add_send ,
                                 const std::vector<EntityKey> & remove_receive,
                                 bool is_full_regen = false);

  bool internal_modification_end( bool regenerate_aura, modification_optimization opt );
  void internal_resolve_shared_modify_delete();
  void internal_resolve_shared_modify_delete_second_pass();
  void internal_resolve_ghosted_modify_delete();
  void internal_resolve_parallel_create();
  void internal_resolve_shared_membership();
  void internal_update_distributed_index( std::vector<Entity> & shared_new );

  /** \brief  Regenerate the shared-entity aura,
   *          adding and removing ghosted entities as necessary.
   *
   *  - a collective parallel operation.
   */
  void internal_regenerate_shared_aura();

  void internal_basic_part_check(const Part* part,
                                 const unsigned ent_rank,
                                 const unsigned undef_rank,
                                 bool& intersection_ok,
                                 bool& rel_target_ok,
                                 bool& rank_ok) const;

  inline void internal_check_unpopulated_relations(Entity entity, EntityRank rank) const;

  // Returns false if there is a problem. It is expected that
  // verify_change_parts will be called if quick_verify_change_part detects
  // a problem, therefore we leave the generation of an exception to
  // verify_change_parts. We want this function to be as fast as
  // possible.
  bool internal_quick_verify_change_part(const Part* part,
                                         const unsigned ent_rank,
                                         const unsigned undef_rank) const;

  void internal_verify_change_parts( const MetaData   & meta ,
                                     const Entity entity ,
                                     const std::vector<Part*> & parts ) const;

  void internal_change_owner_in_comm_data(const EntityKey& key, int new_owner);

  void internal_sync_comm_list_owners();

  void internal_update_fast_comm_maps();

  //------------------------------------

  /** \name  Invariants/preconditions for MetaData.
   * \{
   */

  /** \brief  All non-const methods assert this */
  void require_ok_to_modify() const ;

  void require_entity_owner( const Entity entity, int owner) const ;

  void require_metadata_committed();

  void require_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  bool is_good_rank_and_id(EntityRank ent_rank, EntityId ent_id) const;

  bool is_valid_connectivity(Entity entity, EntityRank rank) const
  {
    if (!is_valid(entity)) return false;
    if (bucket_ptr(entity) == NULL) return false;
    internal_check_unpopulated_relations(entity, rank);
    return true;
  }

  void require_valid_relation( const char action[] ,
                               const BulkData & mesh ,
                               const Entity e_from ,
                               const Entity e_to );


  /** \} */

  //------------------------------------

  // FIXME: Remove this friend once unit-testing has been refactored
  friend class UnitTestModificationEndWrapper;
  friend class ::stk::mesh::MetaData;
  friend class ::stk::mesh::impl::Partition;
  friend class ::stk::mesh::impl::BucketRepository;
  friend class stk::mesh::Bucket; // for field callbacks

#endif /* DOXYGEN_COMPILE */
};

///////////////////////////////////////////////////////////////////////////////
// The following get_connectivity API is designed for *internal* STK_Mesh usage
// to support algorithms that must always be able to get connectivity even
// when it is disabled in the connectivity map. These functions are designed
// to support connectivity-retrieval callsites that need to work regardless of
// the connectivity map. The scratch vectors will be used for allocation when it
// is needed; otherwise, they are ignored. The overloads are provided for the
// common cases where not all connectivity data is needed.
///////////////////////////////////////////////////////////////////////////////

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity, EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinal_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<Permutation> & permutation_scratch_storage );

size_t get_connectivity( const BulkData & mesh,
                         Entity entity,
                         EntityRank to_rank,
                         EntityVector & entity_scratch_storage,
                         std::vector<ConnectivityOrdinal> & ordinal_scratch_storage,
                         std::vector<Permutation> & permutation_scratch_storage );


/** \brief  Is in owned closure of the given process,
 *          typically the local process.
 */
inline
bool in_owned_closure(const BulkData& mesh, const Entity entity , int proc )
{
  const bool same_proc = mesh.parallel_rank() == proc;
  return same_proc && mesh.owned_closure(entity);
}

 /** \brief  Comparitor functor for entities compares the entities' keys */
struct EntityLess {
  EntityLess(const BulkData& mesh) : m_mesh(&mesh) {}

  /** \brief  Comparison operator */
  bool operator()(const Entity lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
    const EntityKey rhs_key = m_mesh->in_index_range(rhs) ? m_mesh->entity_key(rhs) : EntityKey();
    return lhs_key < rhs_key;
  }

  bool operator()(const Entity lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs) ? m_mesh->entity_key(lhs) : EntityKey();
    return lhs_key < rhs;
  }

  bool operator()( const EntityProc & lhs, const EntityProc & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey() ;
    const EntityKey rhs_key = m_mesh->in_index_range(rhs.first) ? m_mesh->entity_key(rhs.first) : EntityKey() ;
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
  }

  bool operator()( const EntityProc & lhs, const Entity rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
    const EntityKey rhs_key = m_mesh->in_index_range(rhs)       ? m_mesh->entity_key(rhs)       : EntityKey();
    return lhs_key < rhs_key;
  }

  bool operator()( const EntityProc & lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = m_mesh->in_index_range(lhs.first) ? m_mesh->entity_key(lhs.first) : EntityKey();
    return lhs_key < rhs ;
  }

  EntityLess& operator=(const EntityLess& rhs)
  {
    m_mesh = rhs.m_mesh;
    return *this;
  }

  const BulkData* m_mesh;
}; //struct EntityLess


inline
BulkData & BulkData::get( const Bucket & bucket) {
  return bucket.bulk_data();
}

inline
BulkData & BulkData::get( const Ghosting & ghost) {
  return ghost.bulk_data();
}

inline
BulkData & BulkData::get( const impl::BucketRepository & bucket_repo ) {
  return bucket_repo.mesh();
}

inline
unsigned BulkData::num_connectivity(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
}

inline
unsigned BulkData::find_ordinal(Entity entity, EntityRank rank, ConnectivityOrdinal ordinal) const
{
  ThrowAssert(bucket_ptr(entity));
  const MeshIndex &mesh_idx = mesh_index(entity);
  unsigned num_rels = mesh_idx.bucket->num_connectivity(mesh_idx.bucket_ordinal, rank);
  ConnectivityOrdinal const *ords = mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);

  unsigned i = 0;
  for (; i < num_rels; ++i)
  {
    if (ords[i] == ordinal)
      break;
  }
  return i;
}

inline
Entity const* BulkData::begin(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::begin_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::begin_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::begin_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::begin_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::begin_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::begin_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->begin_element_permutations(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_nodes(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_edges(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_faces(mesh_idx.bucket_ordinal);
}

inline
unsigned BulkData::num_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->num_elements(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end(mesh_idx.bucket_ordinal, rank);
}

inline
Entity const* BulkData::end_nodes(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_nodes(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_edges(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edges(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_faces(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_faces(mesh_idx.bucket_ordinal);
}

inline
Entity const* BulkData::end_elements(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_elements(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_ordinals(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_ordinals(mesh_idx.bucket_ordinal, rank);
}

inline
ConnectivityOrdinal const* BulkData::end_node_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_edge_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_face_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_ordinals(mesh_idx.bucket_ordinal);
}

inline
ConnectivityOrdinal const* BulkData::end_element_ordinals(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_ordinals(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_permutations(Entity entity, EntityRank rank) const
{
  ThrowAssert(is_valid_connectivity(entity, rank));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_permutations(mesh_idx.bucket_ordinal, rank);
}

inline
Permutation const* BulkData::end_node_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::NODE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_node_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_edge_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::EDGE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_edge_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_face_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::FACE_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_face_permutations(mesh_idx.bucket_ordinal);
}

inline
Permutation const* BulkData::end_element_permutations(Entity entity) const
{
  ThrowAssert(is_valid_connectivity(entity, stk::topology::ELEMENT_RANK));
  const MeshIndex &mesh_idx = mesh_index(entity);
  return mesh_idx.bucket->end_element_permutations(mesh_idx.bucket_ordinal);
}

inline
bool BulkData::has_permutation(Entity entity, EntityRank rank) const
{
  ThrowAssert(bucket_ptr(entity));
  return bucket(entity).has_permutation(rank);
}

/** \} */

inline
void BulkData::internal_basic_part_check(const Part* part,
                                         const unsigned ent_rank,
                                         const unsigned undef_rank,
                                         bool& intersection_ok,
                                         bool& rel_target_ok,
                                         bool& rank_ok) const
{
  // const unsigned part_rank = part->primary_entity_rank();

  intersection_ok = true;
  rel_target_ok   = true;

  // Do we allow arbitrary part changes to entities regardless of part rank? For the sake of the migration, we will for now.
#ifdef SIERRA_MIGRATION
  rank_ok = true;
#else
  const unsigned part_rank = part->primary_entity_rank();
  rank_ok         = ( ent_rank == part_rank ||
                      undef_rank  == part_rank );
#endif
}

inline bool BulkData::internal_quick_verify_change_part(const Part* part,
                                                        const unsigned ent_rank,
                                                        const unsigned undef_rank) const
{
  bool intersection_ok, rel_target_ok, rank_ok;
  internal_basic_part_check(part, ent_rank, undef_rank, intersection_ok, rel_target_ok, rank_ok);
  return intersection_ok && rel_target_ok && rank_ok;
}

inline
int BulkData::entity_comm_owner(const EntityKey & key) const
{
  const int owner_rank = m_entity_comm_map.owner_rank(key);
  ThrowAssertMsg(owner_rank == InvalidProcessRank || owner_rank == parallel_owner_rank(get_entity(key)),
                 "Expect entity " << key.id() << " to have owner " <<
                 parallel_owner_rank(get_entity(key)) << " but in comm map, found " << owner_rank);
  return owner_rank;
}

inline
bool BulkData::in_receive_ghost( EntityKey key ) const
{
  // Ghost communication with owner.
  const int owner_rank = entity_comm_owner(key);
  PairIterEntityComm ec = entity_comm(key);
  return !ec.empty() && ec.front().ghost_id != 0 &&
         ec.front().proc == owner_rank;
}

inline
bool BulkData::in_receive_ghost( const Ghosting & ghost , EntityKey key ) const
{
  const int owner_rank = entity_comm_owner(key);
  return in_ghost( ghost , key , owner_rank );
}

inline
bool BulkData::in_send_ghost( EntityKey key) const
{
  // Ghost communication with non-owner.
  const int owner_rank = entity_comm_owner(key);
  PairIterEntityComm ec = entity_comm(key);
  return ! ec.empty() && ec.back().ghost_id != 0 &&
    ec.back().proc != owner_rank;
}

class LessRelation
{
public:

  LessRelation(const BulkData &mesh) : m_mesh(mesh) { }

  inline bool operator() ( const Relation & lhs , const Relation & rhs ) const
  {
    bool result = false;

    // In Sierra, relations are sorted by RelationType in addition to Rank, Identifier, and target entity key.
  #ifdef SIERRA_MIGRATION
    if (lhs.entity_rank() != rhs.entity_rank()) {
      result = lhs.entity_rank() < rhs.entity_rank();
    }
    else if (lhs.getRelationType() != rhs.getRelationType()) {
      result = lhs.getRelationType() < rhs.getRelationType();
    }
    else if (lhs.relation_ordinal() != rhs.relation_ordinal()) {
      result = lhs.relation_ordinal() < rhs.relation_ordinal();
    }
  #else
    if ( lhs.m_raw_relation.value != rhs.m_raw_relation.value ) {
      result = lhs.m_raw_relation.value < rhs.m_raw_relation.value ;
    }
  #endif
    else {
      Entity lhs_entity = lhs.entity();
      const size_t lhs_offset = m_mesh.is_valid(lhs_entity) ? lhs_entity.local_offset() : Entity::MaxEntity;
      Entity rhs_entity = rhs.entity();
      const size_t rhs_offset = m_mesh.is_valid(rhs_entity) ? rhs_entity.local_offset() : Entity::MaxEntity;
      result = lhs_offset < rhs_offset;
    }
    return result ;
  }

  bool operator() ( const Relation & lhs , Relation::raw_relation_id_type rhs ) const
    { return lhs.raw_relation_id() < rhs ; }

private:

  const BulkData &m_mesh;

  LessRelation();
};

inline
Relation::Relation(const BulkData &mesh,  Entity ent , RelationIdentifier id )
  : m_raw_relation( Relation::raw_relation_id( mesh.entity_rank(ent) , id ) ),
    m_target_entity(ent)
{
#ifdef SIERRA_MIGRATION
  setRelationType(RelationType::INVALID);
#endif
}

inline
bool Bucket::other_entities_have_single_rank(size_type bucket_ordinal, EntityRank rank) const
{
  Entity const * const other_rels = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const other_rels_end = m_dynamic_other_connectivity.end(bucket_ordinal);;

  if ((other_rels == other_rels_end) && (rank != InvalidEntityRank))
    return false;

  return (m_mesh.entity_rank(*other_rels) == rank) && (m_mesh.entity_rank(*(other_rels_end - 1)) == rank);
}


inline
void BulkData::internal_check_unpopulated_relations(Entity entity, EntityRank rank) const
{
#ifndef NDEBUG
  if (m_check_invalid_rels) {
    const MeshIndex &mesh_idx = mesh_index(entity);
    const Bucket &b = *mesh_idx.bucket;
    Bucket::size_type bucket_ord = mesh_idx.bucket_ordinal;
    ThrowAssert(count_valid_connectivity(entity, rank) == b.num_connectivity(bucket_ord, rank));
  }
#endif
}

struct EntityGhostData
{
    enum DIRECTION {
        INVALID,
        NONE,
        SEND,
        RECEIVE
    };
    enum GHOST_LEVEL {
        LOCALLY_OWNED = -1,
        SHARED = 0,
        AURA = 1
    };
    DIRECTION direction;
    int ghostingLevel;
    int processor;
    Entity entity;
    const BulkData * bulkData;

    EntityGhostData()
    : direction(INVALID)
    , ghostingLevel(-2)
    , processor(-1)
    , entity()
    , bulkData(NULL) { }

    // insert to any ostream-like s
    template<class OStream> friend inline OStream& operator << (OStream& s, const DIRECTION& dir)
    {
        switch (dir) {
            case INVALID:
                s << "INVALID";
                break;
            case NONE:
                s << "NONE";
                break;
            case SEND:
                s << "SEND";
                break;
            case RECEIVE:
                s << "RECEIVE";
                break;
            default:
                s << "INVALID";
        }
        return s;
    }
    template<class OStream> inline OStream& printGhostLevel(OStream& s, int gl) const
    {
        switch (gl) {
            case LOCALLY_OWNED:
                s << "LOCALLY_OWNED";
                break;
            case SHARED:
                s << "SHARED";
                break;
            case AURA:
                s << "AURA";
                break;
            default:
                s << "CUSTOM_" << (gl-1);
        }
        return s;
    }
    template<class OStream> friend inline OStream& operator << (OStream& s, const EntityGhostData& egd)
    {
        if (egd.bulkData != NULL) {
            s << "(Entity_gid=";
            s << egd.bulkData->identifier(egd.entity)
              << ", rank=" << egd.bulkData->entity_rank(egd.entity);
        }
        else {
            s << "(Entity_lid=";
            s << egd.entity;
        }
        s << ", direction=" << egd.direction
          << ", processor=" << egd.processor
          << ", ghosting level=";
        egd.printGhostLevel(s,egd.ghostingLevel);
        s << ")";
        return s;
    }
};

void get_ghost_data( const BulkData& bulkData, Entity entity, std::vector<EntityGhostData> & dataVector );



} // namespace mesh
} // namespace stk

#endif //  stk_mesh_BulkData_hpp
