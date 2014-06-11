/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <cstdlib>
#include <stdexcept>

#include <stk_mesh/baseImpl/BucketRepository.hpp>
#include <stk_mesh/baseImpl/Partition.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Trace.hpp>

namespace stk {
namespace mesh {
namespace impl {

BucketRepository::BucketRepository(BulkData & mesh,
                                   unsigned entity_rank_count,
                                   EntityRepository & entity_repo,
                                   const ConnectivityMap& connectivity_map)
  : m_mesh(mesh),
    m_buckets(entity_rank_count),
    m_entity_repo(entity_repo),
    m_partitions(entity_rank_count),
    m_need_sync_from_partitions(entity_rank_count, false),
    m_connectivity_map(connectivity_map)
{
  // Nada.
}

BucketRepository::~BucketRepository()
{
  // Destroy buckets, which were *not* allocated by the set.

  typedef tracking_allocator<Partition, PartitionTag> partition_allocator;

  try {

    for ( std::vector<std::vector<Partition *> >::iterator pv_i = m_partitions.begin();
          pv_i != m_partitions.end(); ++pv_i)
    {
      for (std::vector<Partition *>::iterator p_j = pv_i->begin();
           p_j != pv_i->end(); ++p_j)
      {
        Partition * tmp = *p_j;
        tmp->~Partition();
        partition_allocator().deallocate(tmp,1);
      }
      pv_i->clear();
    }
    m_partitions.clear();
    m_buckets.clear();
  } catch(...) {}
}

size_t BucketRepository::total_field_data_footprint(const FieldBase& f, EntityRank rank) const
{
  if (rank > m_partitions.size())
  {
    return 0;
  }

  size_t retval = 0;
  const std::vector<Partition *> &r_partitions = m_partitions[rank];
  size_t num_partitions = r_partitions.size();
  for (size_t i = 0; i < num_partitions; ++i)
  {
    retval += r_partitions[i]->field_data_footprint(f);
  }
  return retval;
}

void BucketRepository::internal_sort_bucket_entities()
{
  for (std::vector<std::vector<Partition *> >::const_iterator
         i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
  {
    const std::vector<Partition *> & pset = *i ;
    for ( std::vector<Partition*>::const_iterator
            ip = pset.begin() ; ip != pset.end() ; ++ip )
    {
      (*ip)->sort();
    }
  }
}

void BucketRepository::optimize_buckets()
{
  for (std::vector<std::vector<Partition *> >::const_iterator
         i = m_partitions.begin() ; i != m_partitions.end() ; ++i  )
  {
    const std::vector<Partition *> & pset = *i ;
    for ( std::vector<Partition*>::const_iterator
            ip = pset.begin() ; ip != pset.end() ; ++ip )
    {
      (*ip)->compress();
    }
  }
}

////
//// Note that in both versions of get_or_create_partition(..) we need to construct a
//// key vector that the particular format so we can use the lower_bound(..) function to
//// lookup the partition.  Because we are using partitions now instead of buckets, it
//// should be possible to do without that vector and instead do the lookup directly from
//// the PartVector or OrdinalVector.
////

Partition *BucketRepository::get_or_create_partition(
  const unsigned arg_entity_rank ,
  const PartVector &parts)
{
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  TraceIf("stk::mesh::impl::BucketRepository::get_or_create_partition", LOG_BUCKET);

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
                  "Entity rank " << arg_entity_rank << " is invalid");

  // Somehow, this can happen.
  ThrowRequireMsg( !m_buckets.empty(),
                   "m_buckets is empty! Did you forget to initialize MetaData before creating BulkData?");

  std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

  const size_t part_count = parts.size();
  std::vector<unsigned> key(2 + part_count) ;

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , partition_count }
  // Thus partition_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key for a bucket in the partition.
  const unsigned max = static_cast<unsigned>(-1);
  key[0] = part_count+1;
  key[ key[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i )
    {
      key[i+1] = parts[i]->mesh_meta_data_ordinal();
    }
  }

  // If the partition is found, the iterator will be right after it, thanks to the
  // trickiness above.
  const std::vector<Partition *>::iterator ik = lower_bound( partitions , &key[0] );
  const bool partition_exists =
    (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , &key[0] );

  if (partition_exists)
  {
    return ik[-1];
  }


  key[key[0]] = 0;

  typedef tracking_allocator<Partition, PartitionTag> partition_allocator;
  Partition *partition = partition_allocator().allocate(1);
  partition = new (partition) Partition(m_mesh, this, arg_entity_rank, key);

  m_need_sync_from_partitions[arg_entity_rank] = true;
  partitions.insert( ik , partition );

  return partition ;
}

Partition *BucketRepository::get_or_create_partition(
  const unsigned arg_entity_rank ,
  const OrdinalVector &parts)
{
  enum { KEY_TMP_BUFFER_SIZE = 64 };

  TraceIf("stk::mesh::impl::BucketRepository::get_or_create_partition", LOG_BUCKET);

  ThrowRequireMsg(MetaData::get(m_mesh).check_rank(arg_entity_rank),
                  "Entity rank " << arg_entity_rank << " is invalid");

  // Somehow, this can happen.
  ThrowRequireMsg( !m_buckets.empty(),
                   "m_buckets is empty! Did you forget to initialize MetaData before creating BulkData?");

  std::vector<Partition *> & partitions = m_partitions[ arg_entity_rank ];

  const size_t part_count = parts.size();
  std::vector<unsigned> key(2 + part_count) ;

  //----------------------------------
  // Key layout:
  // { part_count + 1 , { part_ordinals } , partition_count }
  // Thus partition_count = key[ key[0] ]
  //
  // for upper bound search use the maximum key for a bucket in the partition.
  const unsigned max = static_cast<unsigned>(-1);
  key[0] = part_count+1;
  key[ key[0] ] = max ;

  {
    for ( unsigned i = 0 ; i < part_count ; ++i ) { key[i+1] = parts[i] ; }
  }

  // If the partition is found, the iterator will be right after it, thanks to the
  // trickiness above.
  const std::vector<Partition *>::iterator ik = lower_bound( partitions , &key[0] );
  const bool partition_exists =
    (ik != partitions.begin()) && raw_part_equal( ik[-1]->key() , &key[0] );

  if (partition_exists)
  {
    return ik[-1];
  }

  key[key[0]] = 0;

  typedef tracking_allocator<Partition, PartitionTag> partition_allocator;
  Partition *partition = partition_allocator().allocate(1);
  partition = new (partition) Partition(m_mesh, this, arg_entity_rank, key);

  m_need_sync_from_partitions[arg_entity_rank] = true;
  partitions.insert( ik , partition );

  return partition ;
}

void BucketRepository::internal_modification_end()
{
  sync_from_partitions();

  // What needs to be done depends on the connectivity map.
  for (EntityRank from_rank = stk::topology::NODE_RANK;
        from_rank < m_connectivity_map.m_map.size();
        ++from_rank)
  {
    const std::vector<Bucket *> &buckets = m_buckets[from_rank];
    unsigned num_buckets = buckets.size();
    for (unsigned j = 0; j < num_buckets; ++j)
    {
      ThrowAssert(buckets[j] != NULL);
      Bucket &bucket = *buckets[j];

      // Update the hop-saving connectivity data on this bucket.
      //
      for (EntityRank to_rank = stk::topology::NODE_RANK;
          to_rank < m_connectivity_map.m_map[from_rank].size();
          ++to_rank)
      {
        switch (m_connectivity_map.m_map[from_rank][to_rank])
        {
        case FIXED_CONNECTIVITY:
          switch (to_rank)
          {
          case stk::topology::NODE_RANK:
            bucket.m_fixed_node_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::EDGE_RANK:
            bucket.m_fixed_edge_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::FACE_RANK:
            bucket.m_fixed_face_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::ELEMENT_RANK:
            bucket.m_fixed_element_connectivity.end_modification(&bucket.m_mesh);
            break;
          default:
            break;
          }
          break;
        case DYNAMIC_CONNECTIVITY:
          switch (to_rank)
          {
          case stk::topology::NODE_RANK:
            bucket.m_dynamic_node_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::EDGE_RANK:
            bucket.m_dynamic_edge_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::FACE_RANK:
            bucket.m_dynamic_face_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::ELEMENT_RANK:
            bucket.m_dynamic_element_connectivity.end_modification(&bucket.m_mesh);
            break;
          case stk::topology::INVALID_RANK:
            break;
          default:
            bucket.m_dynamic_other_connectivity.end_modification(&bucket.m_mesh);
            break;
          }
          break;
        case INVALID_CONNECTIVITY_TYPE:
        default:
          break;
        }
      }
    }
  }
}

void BucketRepository::sync_from_partitions()
{
  for (EntityRank rank = stk::mesh::MetaData::NODE_RANK; rank < m_partitions.size(); ++rank)
  {
    sync_from_partitions(rank);
  }
}

namespace {

inline bool is_null(stk::mesh::impl::Partition *p) { return (p ? false : true);}

}

void BucketRepository::sync_from_partitions(EntityRank rank)
{

  typedef tracking_allocator<Partition, PartitionTag> partition_allocator;

  if (!m_need_sync_from_partitions[rank])
  {
    return;
  }

  std::vector<Partition *> &partitions = m_partitions[rank];

  size_t num_partitions = partitions.size();
  size_t num_buckets = 0;
  for (size_t p_i = 0; p_i < num_partitions; ++p_i)
  {
    if (!partitions[p_i]->empty())
    {
      num_buckets += partitions[p_i]->num_buckets();
    }
  }

  m_buckets[rank].resize(num_buckets);

  bool has_hole = false;
  std::vector<Bucket *>::iterator bkts_i = m_buckets[rank].begin();
  for (size_t p_i = 0; p_i < num_partitions; ++p_i)
  {
    Partition &partition = *partitions[p_i];

    if (partition.empty())
    {
      partitions[p_i]->~Partition();
      partition_allocator().deallocate(partitions[p_i],1);
      partitions[p_i] = 0;
      has_hole = true;
      continue;
    }
    size_t num_bkts_in_partition = partition.num_buckets();
    std::copy(partition.begin(), partition.end(), bkts_i);
    bkts_i += num_bkts_in_partition;
  }

  if (has_hole)
  {
    std::vector<Partition *>::iterator new_end;
    new_end = std::remove_if(partitions.begin(), partitions.end(), is_null);
    size_t new_size = new_end - partitions.begin();  // OK because has_hole is true.
    partitions.resize(new_size);
  }

  sync_bucket_ids(rank);

  m_need_sync_from_partitions[rank] = false;
}

Bucket *BucketRepository::allocate_bucket(EntityRank arg_entity_rank,
                                          const std::vector<unsigned> & arg_key,
                                          size_t arg_capacity )
{

  Bucket * new_bucket = bucket_allocator().allocate(1);
  new_bucket = new (new_bucket) Bucket(m_mesh, arg_entity_rank, arg_key, arg_capacity, m_connectivity_map);
  std::vector<Bucket *> &bucket_vec = m_buckets[arg_entity_rank];
  new_bucket->m_bucket_id = bucket_vec.size();
  bucket_vec.push_back(new_bucket);
  m_need_sync_from_partitions[arg_entity_rank] = true;

  return new_bucket;
}

void BucketRepository::deallocate_bucket(Bucket *b)
{
  ThrowAssertMsg((b != NULL) && (b == m_buckets[b->m_entity_rank][b->m_bucket_id]),
           "BucketRepository::deallocate_bucket(.) m_buckets invariant broken.");

  m_buckets[b->m_entity_rank][b->m_bucket_id] = 0; // space will be reclaimed by sync_from_partitions
  m_need_sync_from_partitions[b->m_entity_rank] = true;
  b->~Bucket();
  bucket_allocator().deallocate(b,1);
}

void BucketRepository::sync_bucket_ids(EntityRank entity_rank)
{
  std::vector<Bucket *> &buckets = m_buckets[entity_rank];
  unsigned num_buckets = buckets.size();
  std::vector<unsigned> id_map(num_buckets);

  for (unsigned i = 0; i < num_buckets; ++i)
  {
    ThrowAssertMsg(buckets[i] != NULL,
                   "BucketRepository::sync_bucket_ids() called when m_buckets["
                   << entity_rank << "] is not dense.");
    id_map[i] = buckets[i]->m_bucket_id;
    buckets[i]->m_bucket_id = i;
  }

  m_mesh.reorder_buckets_callback(entity_rank, id_map);
}

std::vector<Partition *> BucketRepository::get_partitions(EntityRank rank)
{
  if (m_mesh.synchronized_state() != BulkData::SYNCHRONIZED)
  {
    std::vector<Partition *>();
  }
  std::vector<Partition *> retval;
  std::vector<Partition *> &bf_vec = m_partitions[rank];
  for (size_t i = 0; i < bf_vec.size(); ++i)
  {
    retval.push_back(bf_vec[i]);
  }
  return retval;
}

} // namespace impl
} // namespace mesh
} // namespace stk
