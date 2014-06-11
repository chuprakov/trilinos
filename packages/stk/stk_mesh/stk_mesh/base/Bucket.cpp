/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stdlib.h>
#include <memory.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FindRestriction.hpp>

//----------------------------------------------------------------------
namespace stk {
namespace mesh {

namespace {

enum IgnoreMe
{
  DUMMY_VALUE = 0
};

// TODO: When we get C++11, use lambdas instead of these functors

struct CheckSizeFunctor
{
  template <typename Connectivity>
  void operator()(Bucket const& bucket, Connectivity const& connectivity, IgnoreMe) const
  { ThrowAssert(bucket.size() == static_cast<size_t>(connectivity.size())); }

  template <EntityRank Rank, ConnectivityType Type>
  static
  IgnoreMe generate_args(Bucket* other_bucket)
  { return DUMMY_VALUE; }
};

struct AddEntityFunctor
{
  template <typename Connectivity>
  void operator()(Bucket&, Connectivity& connectivity, IgnoreMe)
  { connectivity.add_entity(); }

  template <EntityRank Rank, ConnectivityType Type>
  static
  IgnoreMe generate_args(Bucket* other_bucket)
  { return DUMMY_VALUE; }
};

struct RemoveEntityFunctor
{
  template <typename Connectivity>
  void operator()(Bucket&, Connectivity& connectivity, IgnoreMe)
  { connectivity.remove_entity(); }

  template <EntityRank Rank, ConnectivityType Type>
  static
  IgnoreMe generate_args(Bucket* other_bucket)
  { return DUMMY_VALUE; }
};

struct DeclareRelationFunctor
{
  DeclareRelationFunctor(Bucket::size_type bucket_ordinal, Entity to, ConnectivityOrdinal ordinal,
                         Permutation permutation)
    : m_bucket_ordinal(bucket_ordinal),
      m_to(to),
      m_ordinal(ordinal),
      m_permutation(permutation),
      m_modified(false)
  {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity)
  {
    ThrowAssert( (Connectivity::target_rank == static_cast<EntityRank>(stk::topology::INVALID_RANK) &&
                  bucket.mesh().entity_rank(m_to) > static_cast<EntityRank>(stk::topology::ELEMENT_RANK)) ||
                 bucket.mesh().entity_rank(m_to) == Connectivity::target_rank );
    ThrowAssert(!m_modified);
    m_modified = connectivity.add_connectivity(m_bucket_ordinal, m_to, m_ordinal, m_permutation);
  }

  Bucket::size_type m_bucket_ordinal;
  Entity m_to;
  ConnectivityOrdinal m_ordinal;
  Permutation m_permutation;
  bool m_modified;
};

struct DestroyRelationFunctor
{
  DestroyRelationFunctor(Bucket::size_type bucket_ordinal, Entity to, ConnectivityOrdinal ordinal)
    : m_bucket_ordinal(bucket_ordinal),
      m_to(to),
      m_ordinal(ordinal),
      m_modified(false)
  {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity)
  {
    ThrowAssert( (Connectivity::target_rank == static_cast<EntityRank>(stk::topology::INVALID_RANK) &&
                  bucket.mesh().entity_rank(m_to) > static_cast<EntityRank>(stk::topology::ELEMENT_RANK)) ||
                 bucket.mesh().entity_rank(m_to) == Connectivity::target_rank);
    ThrowAssert(!m_modified);
    m_modified = connectivity.remove_connectivity(m_bucket_ordinal, m_to, m_ordinal);
  }

  Bucket::size_type m_bucket_ordinal;
  Entity m_to;
  ConnectivityOrdinal m_ordinal;
  bool m_modified;
};

struct DebugPrintFunctor
{
  DebugPrintFunctor(std::ostream& out, unsigned ordinal = -1u) : m_out(out), m_ordinal(ordinal) {}

  template <typename Connectivity>
  void operator()(Bucket const& bucket, Connectivity const& connectivity, IgnoreMe)
  { this->operator()(bucket, connectivity); }

  template <typename Connectivity>
  void operator()(Bucket const&, Connectivity const& connectivity)
  {
    if (m_ordinal == -1u) {
      connectivity.debug_dump(m_out);
    }
    else {
      connectivity.debug_dump(m_out, m_ordinal);
    }
  }

  template <EntityRank Rank, ConnectivityType Type>
  static
  IgnoreMe generate_args(Bucket* other_bucket)
  { return DUMMY_VALUE; }

  std::ostream& m_out;
  unsigned m_ordinal;
};

template <typename FixedConnectivity>
void setup_connectivity(stk::topology bucket_topology,
                        EntityRank from_rank,
                        EntityRank to_rank,
                        ConnectivityType& conn_type,
                        FixedConnectivity& fixed_conn,
                        const ConnectivityMap& conn_map)
{
  if (bucket_topology != stk::topology::END_TOPOLOGY && bucket_topology.num_sub_topology(to_rank) > 0 && conn_map(from_rank, to_rank) == FIXED_CONNECTIVITY) {
    fixed_conn.set_num_connectivity(bucket_topology.num_sub_topology(to_rank));
    conn_type = FIXED_CONNECTIVITY;
  }
  else if (from_rank > stk::topology::ELEMENT_RANK || to_rank > stk::topology::ELEMENT_RANK || conn_map(from_rank, to_rank) != INVALID_CONNECTIVITY_TYPE) {
    conn_type = DYNAMIC_CONNECTIVITY;
  }
}

} //namespace anonymous

namespace impl {

struct OverwriteEntityFunctor
{
  OverwriteEntityFunctor(Bucket::size_type old_ordinal, Bucket::size_type new_ordinal) : m_old_ordinal(old_ordinal), m_new_ordinal(new_ordinal) {}

  template <typename Connectivity>
  void operator()(Bucket& bucket, Connectivity& connectivity, Connectivity& old_connectivity)
  { old_connectivity.copy_entity(m_old_ordinal, connectivity, m_new_ordinal); }

  template <EntityRank Rank, ConnectivityType Type>
  static
  impl::BucketConnectivity<Rank, Type>& generate_args(Bucket* other_bucket);

  Bucket::size_type m_old_ordinal;
  Bucket::size_type m_new_ordinal;
};

}

//----------------------------------------------------------------------

bool raw_part_equal( const unsigned * lhs , const unsigned * rhs )
{
  bool result = true ;
  {
    const unsigned * const end_lhs = lhs + *lhs ;
    while ( result && end_lhs != lhs ) {
      result = *lhs == *rhs ;
      ++lhs ; ++rhs ;
    }
  }
  return result ;
}

inline
bool bucket_key_less( const unsigned * lhs , const unsigned * rhs )
{
  const unsigned * const last_lhs = lhs + ( *lhs < *rhs ? *lhs : *rhs );
  while ( last_lhs != lhs && *lhs == *rhs ) { ++lhs ; ++rhs ; }
  return *lhs < *rhs ;
}

// The part count and part ordinals are less
bool BucketLess::operator()( const Bucket * lhs_bucket ,
                             const unsigned * rhs ) const
{ return bucket_key_less( lhs_bucket->key() , rhs ); }

bool BucketLess::operator()( const unsigned * lhs ,
                             const Bucket * rhs_bucket ) const
{ return bucket_key_less( lhs , rhs_bucket->key() ); }

//----------------------------------------------------------------------

Bucket::Bucket( BulkData & arg_mesh ,
                EntityRank arg_entity_rank,
                const std::vector<unsigned> & arg_key,
                size_t arg_capacity,
                const ConnectivityMap& connectivity_map
                )
  : m_mesh(arg_mesh)
  , m_entity_rank(arg_entity_rank)
  , m_topology()
  , m_key(arg_key)
  , m_capacity(arg_capacity)
  , m_size(0)
  , m_bucket(NULL)
  , m_bucket_id(static_cast<unsigned>(-1))
  , m_nodes_per_entity(0)
// TODO: Move owner ranks to BulkData
  , m_entities(arg_capacity)
  , m_owner_ranks(arg_capacity)
  , m_partition(NULL)
  , m_node_kind(INVALID_CONNECTIVITY_TYPE)
  , m_edge_kind(INVALID_CONNECTIVITY_TYPE)
  , m_face_kind(INVALID_CONNECTIVITY_TYPE)
  , m_element_kind(INVALID_CONNECTIVITY_TYPE)
  , m_fixed_node_connectivity()
  , m_fixed_edge_connectivity()
  , m_fixed_face_connectivity()
  , m_fixed_element_connectivity()
  , m_dynamic_node_connectivity(arg_entity_rank, &m_mesh)
  , m_dynamic_edge_connectivity(arg_entity_rank, &m_mesh)
  , m_dynamic_face_connectivity(arg_entity_rank, &m_mesh)
  , m_dynamic_element_connectivity(arg_entity_rank, &m_mesh)
  , m_dynamic_other_connectivity(arg_entity_rank, &m_mesh)
  , m_owned(has_superset(*this, m_mesh.mesh_meta_data().locally_owned_part()))
  , m_shared(has_superset(*this, m_mesh.mesh_meta_data().globally_shared_part()))
{
  ThrowAssertMsg(arg_capacity != 0, "Buckets should never have zero capacity");

  m_topology = get_topology( get_cell_topology(*this), m_mesh.mesh_meta_data().spatial_dimension() );

  if (m_topology != stk::topology::END_TOPOLOGY) {
    m_nodes_per_entity = m_topology.num_nodes();
  }

  setup_connectivity(m_topology, arg_entity_rank, stk::topology::NODE_RANK, m_node_kind, m_fixed_node_connectivity, connectivity_map);
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::EDGE_RANK, m_edge_kind, m_fixed_edge_connectivity, connectivity_map);
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::FACE_RANK, m_face_kind, m_fixed_face_connectivity, connectivity_map);
  setup_connectivity(m_topology, arg_entity_rank, stk::topology::ELEMENT_RANK, m_element_kind, m_fixed_element_connectivity, connectivity_map);

  m_parts.reserve(m_key.size());
  supersets(m_parts);
  m_mesh.new_bucket_callback(m_entity_rank, m_parts, m_capacity);
}

Bucket::~Bucket()
{
  m_mesh.destroy_bucket_callback(m_entity_rank, m_bucket_id, m_capacity);
}

bool Bucket::member( const Part & part ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const unsigned ord = part.mesh_meta_data_ordinal();
  const unsigned * const i = std::lower_bound( i_beg , i_end , ord );

  return i_end != i && ord == *i ;
}

bool Bucket::member_all( const PartVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_all = true ;

  for ( ; result_all && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_all = i_end != i && ord == *i ;
  }
  return result_all ;
}

bool Bucket::member_any( const PartVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const PartVector::const_iterator ip_end = parts.end();
        PartVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = (*ip)->mesh_meta_data_ordinal();
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

bool Bucket::member_any( const OrdinalVector & parts ) const
{
  const unsigned * const i_beg = key() + 1 ;
  const unsigned * const i_end = key() + key()[0] ;

  const OrdinalVector::const_iterator ip_end = parts.end();
        OrdinalVector::const_iterator ip     = parts.begin() ;

  bool result_none = true ;

  for ( ; result_none && ip_end != ip ; ++ip ) {
    const unsigned ord = *ip;
    const unsigned * const i = std::lower_bound( i_beg , i_end , ord );
    result_none = i_end == i || ord != *i ;
  }
  return ! result_none ;
}

unsigned char* Bucket::field_data_location(const FieldBase& field) const
{
  return reinterpret_cast<unsigned char*>(mesh().field_data(field, *this, 0));
}

//----------------------------------------------------------------------

bool has_superset( const Bucket & bucket , const PartVector & ps )
{
  const std::pair<const unsigned *, const unsigned *>
    part_ord = bucket.superset_part_ordinals();

  bool result = ! ps.empty();

  for ( PartVector::const_iterator
        i = ps.begin() ; result && i != ps.end() ; ++i ) {

    const unsigned ordinal = (*i)->mesh_meta_data_ordinal();

    const unsigned * iter =
      std::lower_bound( part_ord.first , part_ord.second , ordinal );

    result = iter < part_ord.second && ordinal == *iter ;
  }
  return result ;
}

void Bucket::supersets( PartVector & ps ) const
{
  const MetaData & mesh_meta_data = MetaData::get( *this );

  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = & mesh_meta_data.get_part( * part_ord.first );
  }
}

void Bucket::supersets( OrdinalVector & ps ) const
{
  std::pair<const unsigned *, const unsigned *>
    part_ord = superset_part_ordinals();

  ps.resize( part_ord.second - part_ord.first );

  for ( unsigned i = 0 ;
        part_ord.first < part_ord.second ; ++(part_ord.first) , ++i ) {
    ps[i] = *part_ord.first;
  }
}

//----------------------------------------------------------------------

bool Bucket::assert_correct() const {
  // test equivalent() method
  const Bucket* bucket = this;
  const Bucket * first = first_bucket_in_partition();
  if (!first || ! bucket->in_same_partition(*first) || ! first->in_same_partition(*bucket) )
    return false;

  // other tests...

  return true;
}

bool Bucket::field_data_is_allocated(const FieldBase& field) const
{
    return mesh().field_is_allocated_for_bucket(field, *this);
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Bucket & k )
{
  const MetaData & mesh_meta_data = MetaData::get(k);
  const std::string & entity_rank_name =
    mesh_meta_data.entity_rank_names()[ k.entity_rank() ];

  const PartVector& parts = k.supersets();

  s << "Bucket( " << entity_rank_name << " : " ;
  for ( PartVector::const_iterator i = parts.begin() ; i != parts.end() ; ++i ) {
    s << (*i)->name() << " " ;
  }
  s << ")" ;

  return s ;
}

std::ostream &
print( std::ostream & os , const std::string & indent , const Bucket & bucket )
{
  const MetaData & mesh_meta_data = MetaData::get(bucket);
  const BulkData & mesh = BulkData::get(bucket);
  const std::string & entity_rank_name =
    mesh_meta_data.entity_rank_names()[ bucket.entity_rank() ];

  const std::pair<const unsigned *, const unsigned *>
    part_ids = bucket.superset_part_ordinals();

  os << "Bucket(size = " << bucket.size() << std::endl << indent << "Part intersection {" ;

  for ( const unsigned * i = part_ids.first ; i < part_ids.second ; ++i ) {
    const Part & part = mesh_meta_data.get_part( *i );
    os << " " << part.name();
  }

  os << " }" << std::endl << indent << entity_rank_name << " members {" ;

  for ( unsigned j = 0 ; j < bucket.size() ; ++j ) {
    const EntityId id = mesh.identifier(bucket[j]);
    os << " " << id ;
  }
  os << " }" << std::endl ;

  bucket.debug_dump(os);
  os << std::endl;

  return os ;
}

struct EntityRankLess
{
  inline bool operator()(Entity entity, EntityRank rank) const {
    return m_mesh->entity_rank(entity) < rank;
  }

  inline bool operator()(EntityRank rank, Entity entity) const {
    return rank < m_mesh->entity_rank(entity);
  }

  const BulkData *m_mesh;
};

size_t Bucket::get_others_begin_index(size_type bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  probe = std::lower_bound(probe, ents_end, rank, cmp);

  return probe - ents_begin;
}

size_t Bucket::get_others_end_index(size_type bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  probe = std::upper_bound(probe, ents_end, rank, cmp);

  return probe - ents_begin;
}

size_t Bucket::get_others_index_count(size_type bucket_ordinal, EntityRank rank) const
{
  Entity const * const ents_begin = m_dynamic_other_connectivity.begin(bucket_ordinal);
  Entity const * const ents_end = m_dynamic_other_connectivity.end(bucket_ordinal);

  EntityRankLess cmp = {&m_mesh};
  Entity const *probe = ents_begin;
  Entity const *const saved_lower = std::lower_bound(probe, ents_end, rank, cmp);
  probe = std::upper_bound(probe, ents_end, rank, cmp);

  return probe - saved_lower;
}

//----------------------------------------------------------------------
// Every bucket in the partition points to the first bucket,
// except the first bucket which points to the last bucket.

Bucket * Bucket::last_bucket_in_partition() const
{
  Bucket * last = last_bucket_in_partition_impl();

  ThrowRequireMsg( NULL != last, "Last is NULL");
  ThrowRequireMsg( last->size() != 0, "Last bucket is empty");

  return last;
}

Bucket * Bucket::last_bucket_in_partition_impl() const
{
  bool this_is_first_bucket_in_partition = (bucket_counter() == 0);

  Bucket * last = NULL;

  if (this_is_first_bucket_in_partition) {
    last = m_bucket;
  } else {
    last = m_bucket->m_bucket;
  }

  return last;
}

//----------------------------------------------------------------------

Bucket * Bucket::first_bucket_in_partition() const
{
  return last_bucket_in_partition_impl()->m_bucket;
}

//----------------------------------------------------------------------

void Bucket::set_last_bucket_in_partition( Bucket * last_bucket )
{
  Bucket * last = last_bucket_in_partition_impl();
  Bucket * first = last->m_bucket;
  first->m_bucket = last_bucket;
}

//----------------------------------------------------------------------

void Bucket::set_first_bucket_in_partition( Bucket * first_bucket )
{
  m_bucket = first_bucket;
}

//----------------------------------------------------------------------

void Bucket::initialize_slot(size_type ordinal, Entity entity)
{
  m_entities[ordinal]    = entity;
  m_owner_ranks[ordinal] = 0;
  if (mesh().is_valid(entity)) {
    mesh().set_state(entity, Created);
  }
}

void Bucket::reset_entity_location(Entity entity, size_type to_ordinal)
{
  Bucket & from_bucket = mesh().bucket(entity);
  const Bucket::size_type from_ordinal = mesh().bucket_ordinal(entity);
  const EntityRank owner_rank = mesh().parallel_owner_rank(entity);

  m_entities[to_ordinal]    = entity;
  m_owner_ranks[to_ordinal] = owner_rank;

  mesh().set_mesh_index(entity, this, to_ordinal);

  m_mesh.copy_entity_fields_callback_same_rank(m_entity_rank, 
                                               m_bucket_id, to_ordinal,
                                              from_bucket.m_bucket_id, from_ordinal);
}

void Bucket::add_entity(Entity entity)
{
  ThrowAssert(m_size < m_capacity);
  ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  ThrowAssert(!mesh().is_valid(entity) || mesh().bucket_ptr(entity) == NULL);
  ThrowAssert(!mesh().is_valid(entity) || mesh().entity_rank(entity) == m_entity_rank);

  initialize_slot(m_size, entity);

  if (mesh().is_valid(entity)) {
    mesh().set_mesh_index(entity, this, m_size);
  }

  ++m_size;

  AddEntityFunctor functor;
  modify_all_connectivity(functor);
}

bool Bucket::destroy_relation(Entity e_from, Entity e_to, const RelationIdentifier local_id )
{
  const size_type from_bucket_ordinal = mesh().bucket_ordinal(e_from);
  DestroyRelationFunctor functor(from_bucket_ordinal, e_to, static_cast<ConnectivityOrdinal>(local_id));
  modify_connectivity(functor, m_mesh.entity_rank(e_to));

  if (functor.m_modified && mesh().bucket(e_from).owned() && (mesh().entity_rank(e_from) > mesh().entity_rank(e_to)) ) {
    --mesh().m_closure_count[e_to.local_offset()];
  }

  return functor.m_modified;
}

bool Bucket::declare_relation(size_type bucket_ordinal, Entity e_to, const ConnectivityOrdinal ordinal, Permutation permutation )
{
  DeclareRelationFunctor functor(bucket_ordinal, e_to, ordinal, permutation);
  modify_connectivity(functor, m_mesh.entity_rank(e_to));

  if (functor.m_modified && owned() && (entity_rank() > mesh().entity_rank(e_to)) ) {
    ++mesh().m_closure_count[e_to.local_offset()];
  }

  return functor.m_modified;
}

void Bucket::remove_entity()
{
  ThrowAssert(m_size > 0);

  --m_size;
  initialize_slot(m_size, Entity());

  RemoveEntityFunctor functor;
  modify_all_connectivity(functor);
}

void Bucket::copy_entity(Entity entity)
{
  ThrowAssert(m_size < m_capacity);
  ThrowAssert(mesh().is_valid(entity));
  ThrowAssert(!mesh().is_valid(m_entities[m_size]));
  ThrowAssert(mesh().bucket_ptr(entity) != NULL);
  ThrowAssert(mesh().bucket_ptr(entity) != this);
  ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  Bucket* old_bucket = mesh().bucket_ptr(entity);
  const Bucket::size_type old_ordinal = mesh().bucket_ordinal(entity);
  reset_entity_location(entity, m_size);

  ++m_size;

  // Unfortunately, we had to copy/paste modify_connectivity to allow dynamic->fixed moves. The
  // modify_connectivity framework couldn't elegantly handle this case.
  switch(m_node_kind) {
  case FIXED_CONNECTIVITY:
    if (old_bucket->m_node_kind == FIXED_CONNECTIVITY) {
      old_bucket->m_fixed_node_connectivity.copy_entity(old_ordinal, m_fixed_node_connectivity);
    }
    else {
      ThrowAssert(old_bucket->m_node_kind != INVALID_CONNECTIVITY_TYPE);
      old_bucket->m_dynamic_node_connectivity.copy_to_fixed(old_ordinal, m_fixed_node_connectivity);
    }
    break;
  case DYNAMIC_CONNECTIVITY: old_bucket->m_dynamic_node_connectivity.copy_entity(old_ordinal, m_dynamic_node_connectivity); break;
  default: break;
  }

  switch(m_edge_kind) {
  case FIXED_CONNECTIVITY:
    if (old_bucket->m_edge_kind == FIXED_CONNECTIVITY) {
      old_bucket->m_fixed_edge_connectivity.copy_entity(old_ordinal, m_fixed_edge_connectivity);
    }
    else {
      ThrowAssert(old_bucket->m_edge_kind != INVALID_CONNECTIVITY_TYPE);
      old_bucket->m_dynamic_edge_connectivity.copy_to_fixed(old_ordinal, m_fixed_edge_connectivity);
    }
    break;
  case DYNAMIC_CONNECTIVITY: old_bucket->m_dynamic_edge_connectivity.copy_entity(old_ordinal, m_dynamic_edge_connectivity); break;
  default: break;
  }

  switch(m_face_kind) {
  case FIXED_CONNECTIVITY:
    if (old_bucket->m_face_kind == FIXED_CONNECTIVITY) {
      old_bucket->m_fixed_face_connectivity.copy_entity(old_ordinal, m_fixed_face_connectivity);
    }
    else {
      ThrowAssert(old_bucket->m_face_kind != INVALID_CONNECTIVITY_TYPE);
      old_bucket->m_dynamic_face_connectivity.copy_to_fixed(old_ordinal, m_fixed_face_connectivity);
    }
    break;
  case DYNAMIC_CONNECTIVITY: old_bucket->m_dynamic_face_connectivity.copy_entity(old_ordinal, m_dynamic_face_connectivity); break;
  default: break;
  }

  switch(m_element_kind) {
  case FIXED_CONNECTIVITY:
    if (old_bucket->m_element_kind == FIXED_CONNECTIVITY) {
      old_bucket->m_fixed_element_connectivity.copy_entity(old_ordinal, m_fixed_element_connectivity);
    }
    else {
      ThrowAssert(old_bucket->m_element_kind != INVALID_CONNECTIVITY_TYPE);
      old_bucket->m_dynamic_element_connectivity.copy_to_fixed(old_ordinal, m_fixed_element_connectivity);
    }
    break;
  case DYNAMIC_CONNECTIVITY: old_bucket->m_dynamic_element_connectivity.copy_entity(old_ordinal, m_dynamic_element_connectivity); break;
  default: break;
  }

  old_bucket->m_dynamic_other_connectivity.copy_entity(old_ordinal, m_dynamic_other_connectivity);
}

void Bucket::overwrite_entity(size_type to_ordinal, Entity entity)
{
  ThrowAssert(to_ordinal < m_capacity);
  ThrowAssert(mesh().is_valid(entity));
  ThrowAssert(mesh().bucket_ptr(entity) != NULL);
  ThrowAssert(mesh().entity_rank(entity) == m_entity_rank);

  const MeshIndex from_index = m_mesh.mesh_index(entity);
  reset_entity_location(entity, to_ordinal);

  impl::OverwriteEntityFunctor functor(from_index.bucket_ordinal, to_ordinal);
  modify_all_connectivity(functor, from_index.bucket);
}

void Bucket::parent_topology( EntityRank parent_rank, std::vector<stk::topology> & parent_topologies) const
{
  PartVector parts;
  supersets(parts);

  parent_topologies.clear();

  for (size_t i=0, e=parts.size(); i<e; ++i)
  {
    Part const & p = *parts[i];
    if ( (p.primary_entity_rank() == parent_rank) && p.topology().is_valid()) {
      parent_topologies.push_back(p.topology());
    }
  }

  std::sort(parent_topologies.begin(),parent_topologies.end());
  std::vector<stk::topology>::iterator end = std::unique(parent_topologies.begin(), parent_topologies.end());

  parent_topologies.erase(end,parent_topologies.end());
}

void Bucket::check_size_invariant() const
{
#ifndef NDEBUG
//  for (size_t i = 0; i < m_entities.size(); ++i) {
//    if (i < m_size) {
//      ThrowAssert(mesh().is_valid(m_entities[i]));
//    }
//    else {
//      ThrowAssert(!mesh().is_valid(m_entities[i]));
//    }
//  }

  CheckSizeFunctor functor;
  const_cast<Bucket*>(this)->modify_all_connectivity(functor);
#endif
}

void Bucket::debug_dump(std::ostream& out, unsigned ordinal) const
{
  DebugPrintFunctor functor(out, ordinal);
  const_cast<Bucket*>(this)->modify_all_connectivity(functor);
}

#ifndef NDEBUG
void Bucket::check_for_invalid_connectivity_request(ConnectivityType const* type) const
{
  EntityRank rank = -1u;
  if (type == &m_node_kind) {
    rank = stk::topology::NODE_RANK;
  }
  else if (type == &m_edge_kind) {
    rank = stk::topology::EDGE_RANK;
  }
  else if (type == &m_face_kind) {
    rank = stk::topology::FACE_RANK;
  }
  else if (type == &m_element_kind) {
    rank = stk::topology::ELEMENT_RANK;
  }
  else {
    ThrowAssert(false);
  }
  // Asking for connectivity between entities of equal rank is always invalid and ok to ask for
  // Asking for connectivity between for FACE_RANK in 2d is always invalid and ok to ask for
  bool isThisEntityAskingForConnectivityToItsOwnRank = entity_rank() == rank;
  bool isThisEntityAskingForFaceConnectivityOnTwoDimensionalMesh = rank == stk::topology::FACE_RANK && mesh().mesh_meta_data().spatial_dimension() == 2;
  ThrowAssert( isThisEntityAskingForConnectivityToItsOwnRank || isThisEntityAskingForFaceConnectivityOnTwoDimensionalMesh);
}
#endif

namespace impl {

template <>
impl::BucketConnectivity<stk::topology::NODE_RANK, FIXED_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::NODE_RANK, FIXED_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_fixed_node_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::EDGE_RANK, FIXED_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::EDGE_RANK, FIXED_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_fixed_edge_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::FACE_RANK, FIXED_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::FACE_RANK, FIXED_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_fixed_face_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::ELEMENT_RANK, FIXED_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::ELEMENT_RANK, FIXED_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_fixed_element_connectivity; }


template <>
impl::BucketConnectivity<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::NODE_RANK, DYNAMIC_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_dynamic_node_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::EDGE_RANK, DYNAMIC_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::EDGE_RANK, DYNAMIC_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_dynamic_edge_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::FACE_RANK, DYNAMIC_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::FACE_RANK, DYNAMIC_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_dynamic_face_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::ELEMENT_RANK, DYNAMIC_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::ELEMENT_RANK, DYNAMIC_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_dynamic_element_connectivity; }

template <>
impl::BucketConnectivity<stk::topology::INVALID_RANK, DYNAMIC_CONNECTIVITY>& OverwriteEntityFunctor::generate_args<stk::topology::INVALID_RANK, DYNAMIC_CONNECTIVITY>(Bucket* other_bucket)
{ return other_bucket->m_dynamic_other_connectivity; }

}

} // namespace mesh
} // namespace stk
