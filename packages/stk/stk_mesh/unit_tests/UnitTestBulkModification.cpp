/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <Shards_BasicTopologies.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_mesh/base/MetaData.hpp>

#include <stk_mesh/fixtures/RingFixture.hpp>

#include <algorithm>
#include <stdexcept>

using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Bucket;
using stk::mesh::BucketIterator;
using stk::mesh::Entity;
using stk::mesh::EntityRank;
using stk::mesh::BucketVector;
using stk::mesh::fixtures::RingFixture;

class UnitTestStkMeshBulkModification {
 public:
  UnitTestStkMeshBulkModification(stk::ParallelMachine pm) :
    m_comm(pm),
    m_num_procs(stk::parallel_machine_size( m_comm )),
    m_rank(stk::parallel_machine_rank( m_comm )),
    m_ring_mesh(pm)
  { }

  void test_bulkdata_not_syncronized();
  void test_closure_of_non_locally_used_entities();
  void test_all_local_nodes();
  void test_all_local_elements();
  void test_parallel_consistency();

  BulkData& initialize_ring_fixture()
  {
    m_ring_mesh.m_meta_data.commit();
    BulkData& bulk_data = m_ring_mesh.m_bulk_data;

    bulk_data.modification_begin();
    m_ring_mesh.generate_mesh( );
    ThrowRequire(bulk_data.modification_end());

    bulk_data.modification_begin();
    m_ring_mesh.fixup_node_ownership( );
    ThrowRequire(bulk_data.modification_end());

    return bulk_data;
  }

  stk::ParallelMachine m_comm;
  int m_num_procs;
  int m_rank;
  RingFixture m_ring_mesh;
};

namespace {

const EntityRank NODE_RANK = stk::mesh::MetaData::NODE_RANK;

STKUNIT_UNIT_TEST( UnitTestBulkDataNotSyrncronized , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_bulkdata_not_syncronized();
}

STKUNIT_UNIT_TEST( UnitTestClosureOfNonLocallyUsedEntities , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_closure_of_non_locally_used_entities();
}

STKUNIT_UNIT_TEST( UnitTestAllLocalNodes , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_all_local_nodes();
}

STKUNIT_UNIT_TEST( UnitTestAllLocalElements , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_all_local_elements();
}

STKUNIT_UNIT_TEST( UnitTestParallelConsistency , testUnit )
{
  UnitTestStkMeshBulkModification unit(MPI_COMM_WORLD);
  unit.test_parallel_consistency();
}

} //end namespace

void UnitTestStkMeshBulkModification::test_bulkdata_not_syncronized()
{
  BulkData& bulk_data = initialize_ring_fixture();

  bulk_data.modification_begin(); // Intentially make things unsynced

  std::vector< Entity> entities;
  std::vector< Entity> entities_closure;
  STKUNIT_ASSERT_THROW(stk::mesh::find_closure(bulk_data, entities, entities_closure), std::runtime_error);
}

void UnitTestStkMeshBulkModification::test_closure_of_non_locally_used_entities()
{
  BulkData& bulk_data = initialize_ring_fixture();

  const stk::mesh::Ghosting & ghost = bulk_data.shared_aura();

  std::vector<stk::mesh::EntityKey> ghost_receive ;

  ghost.receive_list( ghost_receive );

  if (!ghost_receive.empty()) {
    std::vector< Entity> entities;
    std::vector< Entity> entities_closure;

    entities.push_back(bulk_data.get_entity(ghost_receive.front()));

    STKUNIT_ASSERT_THROW(stk::mesh::find_closure(bulk_data, entities, entities_closure), std::runtime_error);
  }
}

void UnitTestStkMeshBulkModification::test_all_local_nodes()
{
  BulkData& bulk_data = initialize_ring_fixture();

  {
    std::vector< Entity> entities;
    std::vector< Entity> entities_closure;
    find_closure(bulk_data, entities, entities_closure);

    // the closure of the an empty set of entities on all procs should be empty
    STKUNIT_EXPECT_TRUE(entities_closure.empty());
  }

  {
    // Get a selector for the univeral part (contains local, shared, and ghosted)
    const stk::mesh::Part& universal = m_ring_mesh.m_meta_data.universal_part();
    stk::mesh::Selector universal_selector(universal);

    // Get the buckets that will give us the universal nodes
    BucketVector buckets = bulk_data.get_buckets(NODE_RANK, universal_selector);

    // Get the universal nodes
    std::vector< Entity> universal_entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    // sort and unique the universal nodes
    std::sort(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityLess(bulk_data));
    std::vector<Entity>::iterator new_end = std::unique(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityEqual());
    universal_entities.erase(new_end, universal_entities.end());

    // Get the buckets that will give us the locally used nodes
    stk::mesh::Selector locally_used_selector =
      m_ring_mesh.m_meta_data.locally_owned_part() |
      m_ring_mesh.m_meta_data.globally_shared_part();

    buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, locally_used_selector);

    // Get the locally used nodes
    std::vector< Entity> entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }

    // Get the closure, passing in the locally used nodes on each proc
    std::vector< Entity> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted nodes on this part will be locally used on one of the other
    // procs, so we expect that they will be part of the closure. In other
    // words, the set of nodes returned by find_closure should exactly match
    // the set of universal nodes.
    STKUNIT_ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    stk::mesh::EntityEqual ee;
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      STKUNIT_EXPECT_TRUE(ee(universal_entities[i], entities_closure[i]));
    }
  }
}

void UnitTestStkMeshBulkModification::test_all_local_elements()
{
  BulkData& bulk_data = initialize_ring_fixture();
  {
    const stk::mesh::Part& universal = m_ring_mesh.m_meta_data.universal_part();
    stk::mesh::Selector universal_selector(universal);

    BucketVector buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, universal_selector);

    // get all the nodes that this process knows about
    std::vector< Entity> universal_entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, universal_selector);

    // get all the elements that this process knows about
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        universal_entities.push_back(*bitr);
      }
    }
    buckets.clear();

    // universal entities should now have all the universal nodes and elements
    // sort and uniq the universal nodes/elements
    std::sort(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityLess(bulk_data));
    std::vector<Entity>::iterator new_end = std::unique(universal_entities.begin(), universal_entities.end(), stk::mesh::EntityEqual());
    universal_entities.erase(new_end, universal_entities.end());

    // get the buckets that we need to traverse to get the locally used elements
    stk::mesh::Selector locally_used_selector =
      m_ring_mesh.m_meta_data.locally_owned_part() |
      m_ring_mesh.m_meta_data.globally_shared_part();

    buckets = bulk_data.get_buckets(stk::topology::ELEMENT_RANK, locally_used_selector);

    // get the locally used elements and store them in entities
    std::vector< Entity> entities;
    for (BucketVector::iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }

    // call find_closure, passing in the locally used elements
    std::vector< Entity> entities_closure;
    stk::mesh::find_closure(bulk_data, entities, entities_closure);

    // The ghosted entities on this proc (element or node) should be contained
    // in the closure of the locally-used element on some other proc, so we
    // expect that they will be part of the closure. In other
    // words, the set of entities returned by find_closure should exactly match
    // the set of universal entities (nodes and elements).
    STKUNIT_ASSERT_TRUE(universal_entities.size() == entities_closure.size());
    stk::mesh::EntityEqual ee;
    for (size_t i = 0; i < entities_closure.size(); ++i) {
      STKUNIT_EXPECT_TRUE(ee(universal_entities[i], entities_closure[i]));
    }
  }
}

void UnitTestStkMeshBulkModification::test_parallel_consistency()
{
  BulkData& bulk_data = initialize_ring_fixture();

  stk::CommBroadcast all(bulk_data.parallel(), 0);

  std::vector< Entity> entities;
  std::vector< Entity> entities_closure;

  // For proc 0 only, add locally used nodes to entities, for all other
  // procs, leave entities empty.
  if (m_rank == 0) {
    stk::mesh::Selector locally_used_selector =
      m_ring_mesh.m_meta_data.locally_owned_part() |
      m_ring_mesh.m_meta_data.globally_shared_part();

    BucketVector const& buckets = bulk_data.get_buckets(stk::topology::NODE_RANK, locally_used_selector);

    for (BucketVector::const_iterator itr = buckets.begin();
         itr != buckets.end(); ++itr) {
      Bucket& b = **itr;
      for (BucketIterator bitr = b.begin(); bitr != b.end(); ++bitr) {
        entities.push_back(*bitr);
      }
    }
  }

  // Call find_closure with proc 0 passing in locally-used nodes
  stk::mesh::find_closure(bulk_data, entities, entities_closure);

  // Proc 0 will broadcast the global ids of the nodes it passed to
  // find_closure

  // pack entities for sizing
  for (std::vector<Entity>::const_iterator
          ep = entities.begin() ; ep != entities.end() ; ++ep ) {
    all.send_buffer().pack<stk::mesh::EntityKey>(bulk_data.entity_key(*ep));
  }

  all.allocate_buffer();

  // pack for real
  for (std::vector<Entity>::const_iterator
         ep = entities.begin() ; ep != entities.end() ; ++ep ) {
    all.send_buffer().pack<stk::mesh::EntityKey>(bulk_data.entity_key(*ep));
  }

  all.communicate();

  // clear-out entities and put the nodes that correspond to the keys
  // broadcast by proc 0 into entities.
  entities.clear();
  stk::CommBuffer& buf = all.recv_buffer();
  stk::mesh::EntityKey k ;
  while ( buf.remaining() ) {
    buf.unpack<stk::mesh::EntityKey>(k);
    Entity e = bulk_data.get_entity(k);
    // If a proc is not aware of a key, that means it has no relationship
    // with that entity, so it can ignore it.
    if (bulk_data.is_valid(e)) {
      entities.push_back(e);
    }
  }

  // sort and unique entities
  std::sort(entities.begin(), entities.end(), stk::mesh::EntityLess(bulk_data));
  std::vector<Entity>::iterator new_end = std::unique(entities.begin(), entities.end(), stk::mesh::EntityEqual());
  entities.erase(new_end, entities.end());

  // If any processor had ghosted nodes that were local to proc 0, those
  // nodes should be in the closure because proc 0 passed them in to
  // find_closure.
  STKUNIT_ASSERT_TRUE(entities.size() == entities_closure.size());
  stk::mesh::EntityEqual ee;
  for (size_t i = 0; i < entities_closure.size(); ++i) {
    STKUNIT_EXPECT_TRUE(ee(entities[i], entities_closure[i]));
  }
}
