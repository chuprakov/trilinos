/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t, NULL
#include <iosfwd>                       // for ostringstream, ostream
#include <set>                          // for set, etc
#include <stdexcept>                    // for runtime_error, logic_error
#include <stk_mesh/base/EntityCommDatabase.hpp>  // for pack_entity_info, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <utility>                      // for pair, make_pair
#include <vector>                       // for vector, etc
#include "Shards_BasicTopologies.hpp"   // for getCellTopologyData, etc
#include "gtest/gtest.h"                // for AssertHelper
#include "mpi.h"                        // for MPI_COMM_WORLD, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, Bucket::iterator
#include "stk_mesh/base/BulkData.hpp"   // for BulkData, EntityLess, etc
#include "stk_mesh/base/CellTopology.hpp"  // for CellTopology
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting, operator<<
#include "stk_mesh/base/Types.hpp"      // for EntityProc, PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/Parallel.hpp"  // for parallel_machine_size, etc
#include "stk_util/parallel/ParallelComm.hpp"  // for CommAll, CommBuffer
namespace stk { namespace mesh { class Part; } }

using stk::mesh::Part;
using stk::mesh::MetaData;

// UnitTestBulkData_new is the beginnings of a refactoring of the bulk
// data unit test.  It relies on a customized BoxFixture to rapidly
// create a mesh for testing.

namespace {

void new_insert_transitive_closure( stk::mesh::BulkData& bulk_data, std::set<stk::mesh::EntityProc,stk::mesh::EntityLess> &  ,
					 const stk::mesh::EntityProc & entry );
void new_comm_sync_send_recv(
   stk::mesh::BulkData & mesh ,
   std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > & new_send ,
   std::set< stk::mesh::Entity , stk::mesh::EntityLess > & new_recv );

void new_comm_recv_to_send(
  stk::mesh::BulkData & mesh ,
  const std::set< stk::mesh::Entity , stk::mesh::EntityLess > & new_recv ,
        std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > & new_send );

/**
 * The customized box fixture used in this file for testing. This fixture
 * is similar to the BoxFixture it inherits from, with the only difference
 * being the extra parts that this fixture declares for testing purposes.
 */
struct TestBoxFixture : public stk::mesh::fixtures::BoxFixture
{
  TestBoxFixture(stk::ParallelMachine pm = MPI_COMM_WORLD,
                 unsigned block_size = 1000) :
    BoxFixture(pm, block_size),
    m_test_part ( m_fem_meta.declare_part ( "Test Part" ) ),
    m_cell_part ( m_fem_meta.declare_part ( "Cell list" , stk::topology::ELEM_RANK ) ),
    m_part_A_0 ( m_fem_meta.declare_part ( "Part A 0", stk::topology::NODE_RANK ) ),
    m_part_A_1 ( m_fem_meta.declare_part ( "Part A 1", stk::topology::EDGE_RANK ) ),
    m_part_A_2 ( m_fem_meta.declare_part ( "Part A 2", stk::topology::FACE_RANK ) ),
    m_part_A_3 ( m_fem_meta.declare_part ( "Part A 3", stk::topology::ELEM_RANK ) ),
    m_part_A_superset ( m_fem_meta.declare_part ( "Part A superset" ) ),
    m_part_B_0 ( m_fem_meta.declare_part ( "Part B 0", stk::topology::NODE_RANK ) ),
    m_part_B_1 ( m_fem_meta.declare_part ( "Part B 1", stk::topology::EDGE_RANK ) ),
    m_part_B_2 ( m_fem_meta.declare_part ( "Part B 2", stk::topology::FACE_RANK ) ),
    m_part_B_3 ( m_fem_meta.declare_part ( "Part B 3", stk::topology::ELEM_RANK ) ),
    m_part_B_superset ( m_fem_meta.declare_part ( "Part B superset" ) )
  {
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_0 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_1 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_2 );
    m_fem_meta.declare_part_subset ( m_part_A_superset , m_part_A_3 );

    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_0 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_1 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_2 );
    m_fem_meta.declare_part_subset ( m_part_B_superset , m_part_B_3 );

    // None of the tests currently need to make any addtional changes
    // to MetaData; if this changes, the line below will have to be
    // removed.
    m_fem_meta.commit();
  }

  Part     & m_test_part;   // A simple part
  Part     & m_cell_part;   // A part to put cells in

  Part     & m_part_A_0;
  Part     & m_part_A_1;
  Part     & m_part_A_2;
  Part     & m_part_A_3;

  Part     & m_part_A_superset;

  Part     & m_part_B_0;
  Part     & m_part_B_1;
  Part     & m_part_B_2;
  Part     & m_part_B_3;

  Part     & m_part_B_superset;
};

}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyAssertOwnerDeletedEntity )
{
  TestBoxFixture fixture;

  stk::mesh::BulkData         &bulk = fixture.bulk_data();
  stk::mesh::Part             &new_part = fixture.m_test_part;
  stk::mesh::PartVector        add_part;
  add_part.push_back ( &new_part );

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  // Find a cell owned by this process
  stk::mesh::Entity cell_to_delete = stk::mesh::Entity();
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) == fixture.comm_rank() )
      {
        cell_to_delete = *cur_entity;
        break;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  STKUNIT_ASSERT ( bulk.is_valid(cell_to_delete) );
  bulk.modification_begin();
  bulk.destroy_entity ( cell_to_delete );
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDetectsBadKey )
{
  TestBoxFixture fixture;

  stk::mesh::BulkData         &bulk = fixture.bulk_data();
  stk::mesh::Part             &new_part = fixture.m_test_part;
  stk::mesh::PartVector        add_part, empty_vector;
  add_part.push_back ( &new_part );

  stk::mesh::EntityKey bad_key1 ( static_cast<stk::mesh::EntityRank>(45) , 1 );  // Bad entity rank
  stk::mesh::EntityKey bad_key2 ( stk::topology::EDGE_RANK , 0 );   // Bad id

  STKUNIT_ASSERT_THROW ( bulk.declare_entity(bad_key1.rank(),
                                             bad_key1.id(),
                                             empty_vector),
                         std::logic_error );
  STKUNIT_ASSERT_THROW ( bulk.declare_entity(bad_key2.rank(),
                                             bad_key2.id(),
                                             empty_vector),
                         std::logic_error );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDetectsNonOwnerChange )
{
  // Set up a mesh where there are shared nodes. Take one of the nodes, and
  // have the non-owning processes try to make a change to that node; this
  // should cause an exception.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  int p_size = stk::parallel_machine_size(pm);
  int p_rank = stk::parallel_machine_rank(pm);

  stk::mesh::fixtures::QuadFixture fixture(pm, 1 /*nx*/, p_size /*ny*/);
  fixture.m_meta.commit();
  fixture.generate_mesh();
  stk::mesh::BulkData & bulk = fixture.m_bulk_data;

  stk::mesh::PartVector empty_vector;

  stk::mesh::Entity shared_node = fixture.node(1 /*x*/, 1 /*y*/);
  // Assert that this node is shared
  if ( p_size > 1 && bulk.is_valid(shared_node) && (p_rank == 0 || p_rank == 1) ) {
    STKUNIT_ASSERT_GE(bulk.entity_comm_sharing(bulk.entity_key(shared_node)).size(), 1u);
  }

  bulk.modification_begin();

  // Non-owners of shared_node will attempt to make a change to it; this should
  // cause an exception
  if (bulk.is_valid(shared_node) && p_rank != bulk.parallel_owner_rank(shared_node)) {
    STKUNIT_ASSERT_THROW(bulk.change_entity_parts(shared_node,
                                                  empty_vector,  //add parts
                                                  empty_vector), //rem parts
                         std::logic_error);
  }

  bulk.modification_end();
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyExplicitAddInducedPart )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData     &bulk = fixture.bulk_data ();
  stk::mesh::PartVector    empty_vector;
  stk::mesh::PartVector    cell_part_vector;

  bulk.modification_begin();

  stk::mesh::Entity new_cell = bulk.declare_entity ( stk::topology::ELEMENT_RANK , fixture.comm_rank()+1 , empty_vector );
  stk::mesh::Entity new_node = bulk.declare_entity ( stk::topology::NODE_RANK , fixture.comm_rank()+1 , empty_vector );

  bulk.declare_relation ( new_cell , new_node , 1 );

  cell_part_vector.push_back ( &fixture.m_cell_part );
  bulk.change_entity_parts ( new_cell , cell_part_vector );
#ifdef SIERRA_MIGRATION
  bulk.change_entity_parts ( new_node , cell_part_vector );
#else
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_node , cell_part_vector ) , std::runtime_error );
#endif
}

/************************
 * This unit test is not possible currently because of the lack of
 * separation between internal part modification routines and public
 * part modification routines.
STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCannotRemoveFromSpecialParts )
{
  stk::mesh::fixtures::BoxFixture fixture;
  stk::mesh::BulkData          &bulk = fixture.bulk_data();
  stk::mesh::PartVector         test_parts;
  stk::mesh::PartVector         out_parts;
  stk::mesh::PartVector         empty_vector;

  stk::mesh::Entity new_cell = bulk.declare_entity ( 3 , fixture.comm_rank()+1 , empty_vector );
  test_parts.push_back ( &fixture.fem_meta().universal_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.fem_meta().locally_owned_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
  test_parts.clear();
  test_parts.push_back ( &fixture.fem_meta().globally_shared_part() );
  STKUNIT_ASSERT_THROW ( bulk.change_entity_parts ( new_cell , empty_vector , test_parts ) , std::runtime_error );
}
 */


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyDefaultPartAddition )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData            &bulk = fixture.bulk_data ();

  bulk.modification_begin();
  stk::mesh::Entity new_cell = fixture.get_new_entity ( stk::topology::ELEM_RANK , 1 );
  bulk.modification_end();

  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.fem_meta().universal_part() ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.fem_meta().locally_owned_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangePartsSerial )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData            &bulk = fixture.bulk_data ();
  stk::mesh::PartVector           create_parts , remove_parts , add_parts, empty_parts;

  create_parts.push_back ( &fixture.m_test_part );
  create_parts.push_back ( &fixture.m_part_A_3 );
  remove_parts.push_back ( &fixture.m_part_A_3 );
  add_parts.push_back ( &fixture.m_part_B_superset );
  add_parts.push_back ( &fixture.m_cell_part );

  bulk.modification_begin();
  stk::mesh::Entity new_cell = fixture.get_new_entity ( stk::topology::ELEM_RANK , 1 );
  bulk.change_entity_parts ( new_cell , create_parts , empty_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , add_parts , remove_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.change_entity_parts ( new_cell , empty_parts , add_parts );
  bulk.modification_end();
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_test_part ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_part_A_3 ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.m_part_A_superset ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_part_B_superset ) );
  STKUNIT_ASSERT ( !bulk.bucket(new_cell).member ( fixture.m_cell_part ) );

  //Verify still a member of default parts
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.fem_meta().universal_part() ) );
  STKUNIT_ASSERT ( bulk.bucket(new_cell).member ( fixture.fem_meta().locally_owned_part() ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyParallelAddParts )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData             &bulk = fixture.bulk_data ();
  stk::mesh::PartVector            add_part;

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  add_part.push_back ( &fixture.m_part_A_0 );

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i =  bulk.comm_list().begin();
        i != bulk.comm_list().end() ; ++i ) {
    if ( i->key.rank() == 0 ) {
      if ( i->owner == fixture.comm_rank() ) {
        bulk.change_entity_parts ( i->entity, add_part, stk::mesh::PartVector() );
      }
    }
  }

  bulk.modification_end();

  for ( stk::mesh::EntityCommListInfoVector::const_iterator
        i =  bulk.comm_list().begin();
        i != bulk.comm_list().end() ; ++i ) {
    if ( i->key.rank() == 0 ) {
      STKUNIT_ASSERT ( bulk.bucket(i->entity).member ( fixture.m_part_A_0 ) );
    }
  }
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyInducedMembership )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData             &bulk = fixture.bulk_data ();
  stk::mesh::PartVector            create_node_parts , create_cell_parts , empty_parts;

  create_node_parts.push_back ( &fixture.m_part_A_0 );
  create_cell_parts.push_back ( &fixture.m_cell_part );

  bulk.modification_begin();

  stk::mesh::Entity node = fixture.get_new_entity ( stk::topology::NODE_RANK , 1 );
  stk::mesh::Entity cell = fixture.get_new_entity ( stk::topology::ELEM_RANK , 1 );

  bulk.modification_begin();

  bulk.change_entity_parts ( node , create_node_parts , stk::mesh::PartVector () );
  bulk.change_entity_parts ( cell , create_cell_parts , stk::mesh::PartVector () );
  // Add node to cell part
  stk::mesh::RelationIdentifier cell_node_rel_id = 0;
  bulk.declare_relation ( cell , node , cell_node_rel_id );
  bulk.modification_end();

  STKUNIT_ASSERT ( bulk.bucket(node).member ( fixture.m_cell_part ) );

  bulk.modification_begin();
  bulk.destroy_relation ( cell , node, cell_node_rel_id );
  bulk.modification_end();

  STKUNIT_ASSERT ( !bulk.bucket(node).member ( fixture.m_cell_part ) );
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCanRemoveFromSetWithDifferentRankSubset )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData           &bulk = fixture.bulk_data ();
  stk::mesh::PartVector          add_parts , remove_parts, empty_parts;

  add_parts.push_back ( &fixture.m_part_B_3 );
  add_parts.push_back ( &fixture.m_part_A_superset );

  remove_parts.push_back ( &fixture.m_part_A_superset );

  bulk.modification_begin();

  stk::mesh::Entity e = bulk.declare_entity ( stk::topology::ELEMENT_RANK , fixture.comm_rank()+1 , add_parts );
  bulk.modification_end();

  bulk.modification_begin();
  bulk.change_entity_parts ( e , empty_parts , remove_parts );
  bulk.modification_end();

  STKUNIT_ASSERT ( bulk.bucket(e).member ( fixture.m_part_B_3 ) );
  STKUNIT_ASSERT ( !bulk.bucket(e).member ( fixture.m_part_A_superset ) );
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyCommonGhostingName )
{

  TestBoxFixture fixture;
  stk::mesh::BulkData          &bulk = fixture.bulk_data ();

  bulk.modification_begin();

  if ( fixture.comm_size() == 1 ) return;

  if ( fixture.comm_rank() == 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 1" ) , std::runtime_error );
  }
  else
  {
    STKUNIT_ASSERT_THROW ( bulk.create_ghosting ( "Name 2" ) , std::runtime_error );
  }
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyTrivialDestroyAllGhostings )
{
  TestBoxFixture fixture;

  if ( fixture.comm_size() == 1 ) return;

  stk::mesh::BulkData  &bulk = fixture.bulk_data();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  stk::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );

  // Find a cell owned by this process
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;

  std::vector<stk::mesh::EntityProc>  to_send;
  std::vector<stk::mesh::EntityKey>   empty_vector;
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) == fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        send_rank++;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }
  bulk.change_ghosting ( ghosting , to_send , empty_vector );
  bulk.modification_end();


  {
    std::vector<stk::mesh::EntityProc> send_list ;
    std::vector<stk::mesh::EntityKey>  recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    STKUNIT_ASSERT ( ! send_list.empty()  );
    STKUNIT_ASSERT ( ! recv_list.empty() );
  }

  // Usage of operator << in Ghosting.cpp
  std::ostringstream oss;
  oss << ghosting;

  bulk.modification_begin();
  bulk.destroy_all_ghosting ();
  bulk.modification_end();

  {
    std::vector<stk::mesh::EntityProc> send_list ;
    std::vector<stk::mesh::EntityKey>  recv_list ;
    ghosting.send_list( send_list );
    ghosting.receive_list( recv_list );

    STKUNIT_ASSERT ( send_list.empty() );
    STKUNIT_ASSERT ( recv_list.empty() );
  }
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyChangeGhostingGuards )
{
  TestBoxFixture fixture1, fixture2;
  stk::mesh::BulkData & bulk1 = fixture1.bulk_data ();
  stk::mesh::BulkData & bulk2 = fixture2.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box1[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };
  int local_box2[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk1.modification_begin();
  fixture1.generate_boxes( root_box, local_box1 );
  STKUNIT_ASSERT(bulk1.modification_end());

  bulk2.modification_begin();
  fixture2.generate_boxes( root_box, local_box2 );
  STKUNIT_ASSERT(bulk2.modification_end());

  bulk1.modification_begin();
  bulk2.modification_begin();

  std::vector<stk::mesh::EntityProc>  to_send;
  std::vector<stk::mesh::EntityKey>   empty_vector;
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk1.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;
  while ( cur_bucket != bulk1.buckets(stk::topology::ELEM_RANK).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk1.parallel_owner_rank(*cur_entity) == fixture1.comm_rank() )
      {
        if ( send_rank == fixture1.comm_size() ) send_rank = 0;
        if ( send_rank != fixture1.comm_rank() )
          to_send.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        ++send_rank;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  stk::mesh::Ghosting &ghosting = bulk1.create_ghosting ( "Ghost 1" );
  STKUNIT_ASSERT_THROW ( bulk1.change_ghosting ( bulk1.shared_aura() , to_send , empty_vector ) , std::runtime_error );

  ghosting.receive_list(empty_vector);
  ghosting.send_list(to_send);

  bulk1.modification_end();
  bulk2.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyOtherGhostingGuards )
{
  TestBoxFixture fixture;
  stk::mesh::BulkData          &bulk = fixture.bulk_data ();

  const int root_box[3][2] = { { 0 , 4 } , { 0 , 5 } , { 0 , 6 } };
  int local_box[3][2] = { { 0 , 0 } , { 0 , 0 } , { 0 , 0 } };

  bulk.modification_begin();
  fixture.generate_boxes( root_box, local_box );
  STKUNIT_ASSERT(bulk.modification_end());

  bulk.modification_begin();

  std::vector<stk::mesh::EntityProc>  to_send_unowned;
  std::vector<stk::mesh::EntityProc>  empty_send;
  std::vector<stk::mesh::EntityKey>   to_remove_not_ghosted;
  std::vector<stk::mesh::EntityKey>   empty_remove;
  std::vector<stk::mesh::Bucket *>::const_iterator cur_bucket = bulk.buckets(stk::topology::ELEM_RANK).begin();
  int send_rank = 0;
  while ( cur_bucket != bulk.buckets(stk::topology::ELEM_RANK).end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) != fixture.comm_rank() )
      {
        if ( send_rank == fixture.comm_size() ) send_rank = 0;
        if ( send_rank != fixture.comm_rank() )
          to_send_unowned.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        ++send_rank;
      }
      else
      {
        to_remove_not_ghosted.push_back ( bulk.entity_key(*cur_entity) );
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  stk::mesh::Ghosting &ghosting = bulk.create_ghosting ( "Ghost 1" );
  if ( to_send_unowned.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , to_send_unowned , empty_remove );
  }

  if ( to_remove_not_ghosted.size() > 0 )
  {
    STKUNIT_ASSERT_THROW ( bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted ) , std::runtime_error );
  }
  else
  {
    bulk.change_ghosting ( ghosting , empty_send , to_remove_not_ghosted );
  }
  bulk.modification_end();
}


STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyPartsOnCreate )
{
   TestBoxFixture fixture;
   stk::mesh::BulkData           & bulk = fixture.bulk_data ();
   stk::mesh::Part               & part_a = fixture.m_part_A_0;
   stk::mesh::Part               & part_b = fixture.m_part_B_0;

   stk::mesh::PartVector           create_vector;
   create_vector.push_back ( &part_a );

   bulk.modification_begin();

   stk::mesh::Entity node = bulk.declare_entity ( stk::topology::NODE_RANK , fixture.comm_rank()+1 ,create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( bulk.bucket(node).member ( part_a ) );

   bulk.modification_begin();
   create_vector.push_back ( &part_b );
   stk::mesh::Entity node2 = bulk.declare_entity ( stk::topology::NODE_RANK , fixture.comm_size() + fixture.comm_rank() + 1 , create_vector );
   bulk.modification_end();

   STKUNIT_ASSERT ( bulk.bucket(node2).member ( part_a ) );
   STKUNIT_ASSERT ( bulk.bucket(node2).member ( part_b ) );
}

//----------------------------------------------------------------------

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , verifyBoxGhosting )
{
  const int p_size = stk::parallel_machine_size( MPI_COMM_WORLD );
  if ( 8 < p_size ) { return ; }

  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, 2, 2, 2 );
  fixture.m_meta.commit();
  fixture.generate_mesh();
  const stk::mesh::BulkData& mesh = fixture.m_bulk_data;

  for ( size_t iz = 0 ; iz < 3 ; ++iz ) {
    for ( size_t iy = 0 ; iy < 3 ; ++iy ) {
      for ( size_t ix = 0 ; ix < 3 ; ++ix ) {
        stk::mesh::Entity const node = fixture.node(ix,iy,iz);
        STKUNIT_ASSERT( mesh.is_valid(node) );

        STKUNIT_ASSERT( fixture.node_id(ix,iy,iz) == mesh.identifier(node) );
        stk::mesh::fixtures::HexFixture::Scalar * const node_coord =
            stk::mesh::field_data( fixture.m_coord_field , node );
        STKUNIT_ASSERT( node_coord != NULL );
      }
    }
  }

  for ( size_t iz = 0 ; iz < 2 ; ++iz ) {
  for ( size_t iy = 0 ; iy < 2 ; ++iy ) {
  for ( size_t ix = 0 ; ix < 2 ; ++ix ) {
    stk::mesh::Entity const elem = fixture.elem(ix,iy,iz);
    STKUNIT_ASSERT( mesh.is_valid(elem) );
    size_t num_elem_nodes = mesh.num_nodes(elem);
    STKUNIT_ASSERT_EQUAL( 8u , num_elem_nodes );
    stk::mesh::Entity const *elem_nodes = mesh.begin_nodes(elem);
    // stk::mesh::ConnectivityOrdinal const *elem_node_ords = mesh.begin_node_ordinals(elem);
    if ( 8u == num_elem_nodes ) {
      STKUNIT_ASSERT( elem_nodes[0] == fixture.node(ix,iy,iz));
      STKUNIT_ASSERT( elem_nodes[1] == fixture.node(ix+1,iy,iz));
      STKUNIT_ASSERT( elem_nodes[2] == fixture.node(ix+1,iy+1,iz));
      STKUNIT_ASSERT( elem_nodes[3] == fixture.node(ix,iy+1,iz));
      STKUNIT_ASSERT( elem_nodes[4] == fixture.node(ix,iy,iz+1));
      STKUNIT_ASSERT( elem_nodes[5] == fixture.node(ix+1,iy,iz+1));
      STKUNIT_ASSERT( elem_nodes[6] == fixture.node(ix+1,iy+1,iz+1));
      STKUNIT_ASSERT( elem_nodes[7] == fixture.node(ix,iy+1,iz+1));
    }
    // Now check access to field data via the fast rank functions.
    // stk::mesh::Node const *eph_elem_nodes = mesh.begin_nodes(elem);
    // for ( size_t j = 0 ; j < num_elem_nodes ; ++j )
    // {
    //   stk::mesh::fixtures::HexFixture::Scalar * const node_coord =
    //     stk::mesh::field_data( fixture.m_coord_field , eph_elem_nodes[j]);
    //   STKUNIT_EXPECT_EQ( node_coord, elem_node_coord[ elem_node_ords[j] ] );
    // }

  }
  }
  }
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , testEntityComm )
{
  //Test on unpack_field_values in EntityComm.cpp
  //code based on ../base/BulkDataGhosting.cpp
  //Create a simple mesh. Add nodes one element and some parts.

  const int spatial_dimension = 3;

  stk::mesh::MetaData fem_meta(spatial_dimension);

  stk::mesh::CellTopology tet_top(shards::getCellTopologyData<shards::Tetrahedron<4> >());
  stk::mesh::Part & part_a = fem_meta.declare_part( "block_a", tet_top );
  stk::mesh::Part & part_b = fem_meta.declare_part( "block_b", tet_top );

  stk::mesh::CellTopology node_top(shards::getCellTopologyData<shards::Node>());
  stk::mesh::Part & part_a_0 = fem_meta.declare_part( "block_a_0", node_top );

  typedef stk::mesh::Field<double>  ScalarFieldType;

  ScalarFieldType & volume =
    fem_meta.declare_field < ScalarFieldType > ( stk::topology::ELEMENT_RANK, "volume" , 4 );
  ScalarFieldType & temperature =
    fem_meta.declare_field < ScalarFieldType > ( stk::topology::ELEMENT_RANK, "temperature" , 4 );
  stk::mesh::Part  & universal     = fem_meta.universal_part ();
  put_field ( volume , universal );
  put_field ( temperature , universal );

  fem_meta.commit();

  stk::mesh::PartVector    create_vector;
  stk::mesh::PartVector    empty_vector;
  create_vector.push_back ( &part_a );
  create_vector.push_back ( &part_b );

  stk::mesh::BulkData bulk ( fem_meta , MPI_COMM_WORLD , 100 );

  bulk.modification_begin();

  stk::mesh::Ghosting &ghosts = bulk.create_ghosting ( "Ghost 1" );

  int size2 = stk::parallel_machine_size( MPI_COMM_WORLD );
  int rank_count2 = stk::parallel_machine_rank( MPI_COMM_WORLD );
  int new_id2 = size2 + rank_count2;

  stk::mesh::Entity elem2 = bulk.declare_entity ( stk::topology::ELEMENT_RANK , new_id2+1 ,create_vector );
  STKUNIT_ASSERT_EQUAL( bulk.bucket(elem2).member ( part_a ), true );

  int size = stk::parallel_machine_size( MPI_COMM_WORLD );
  int rank_count = stk::parallel_machine_rank( MPI_COMM_WORLD );

  int id_base = 0;
  for ( id_base = 0 ; id_base < 99 ; ++id_base )
  {
    int new_id = size * id_base + rank_count;
    stk::mesh::Entity new_node = bulk.declare_entity( stk::topology::NODE_RANK , new_id+1 , empty_vector );
    STKUNIT_ASSERT_EQUAL( bulk.bucket(new_node).member ( part_a_0 ), false );
  }

  //Create a bucket of nodes for sending

  std::vector<stk::mesh::EntityProc>  add_send;

  const std::vector<stk::mesh::Bucket*> & buckets = bulk.buckets( stk::topology::NODE_RANK );

  std::vector<stk::mesh::Bucket*>::const_iterator cur_bucket;

  cur_bucket = buckets.begin();

  int send_rank = 0;
  while ( cur_bucket != buckets.end() )
  {
    stk::mesh::Bucket::iterator cur_entity = (*cur_bucket)->begin();
    while ( cur_entity != (*cur_bucket)->end() )
    {
      if ( bulk.parallel_owner_rank(*cur_entity) == rank_count )
      {
        if ( send_rank == size ) send_rank = 0;
        if ( send_rank != rank_count )
          add_send.push_back ( std::make_pair ( *cur_entity , send_rank ) );
        ++send_rank;
      }
      ++cur_entity;
    }
    ++cur_bucket;
  }

  stk::mesh::EntityLess entless(bulk);
  std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > new_send(entless) ;
  std::set< stk::mesh::Entity ,   stk::mesh::EntityLess > new_recv(entless) ;

  //  Keep the closure of the remaining received ghosts.
  //  Working from highest-to-lowest key (rank entity type)
  //  results in insertion of the transitive closure.
  //  Insertion will not invalidate the associative container's iterator.

  for ( std::set< stk::mesh::Entity , stk::mesh::EntityLess >::iterator
        i = new_recv.end() ; i != new_recv.begin() ; ) {
    --i ;

    const unsigned erank = bulk.entity_rank(*i);

    stk::mesh::MeshIndex mesh_idx = bulk.mesh_index(*i);
    stk::mesh::Bucket &bkt = *mesh_idx.bucket;
    stk::mesh::Ordinal bkt_ordinal = mesh_idx.bucket_ordinal;

    for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK;
          irank < erank;
          ++irank)
    {
      stk::mesh::Entity const *irels_itr = bkt.begin(bkt_ordinal, irank);
      stk::mesh::Entity const *irels_end = bkt.end(bkt_ordinal, irank);
      for (; irels_itr != irels_end; ++irels_itr)
      {
        if (bulk.in_receive_ghost( ghosts , bulk.entity_key(*irels_itr) ) )
        {
          new_recv.insert( *irels_itr );
        }
      }
    }
  }

  //  Initialize the new_send from the new_recv
  new_comm_recv_to_send( bulk , new_recv , new_send );

  //------------------------------------
  // Add the specified entities and their closure to the send ghosting

  for ( std::vector< stk::mesh::EntityProc >::const_iterator
        i = add_send.begin() ; i != add_send.end() ; ++i ) {
        new_insert_transitive_closure( bulk, new_send , *i );
  }

  // Synchronize the send and receive list.
  // If the send list contains a not-owned entity
  // inform the owner and receiver to ad that entity
  // to their ghost send and receive lists.

  new_comm_sync_send_recv( bulk , new_send , new_recv );

  //------------------------------------
  // Push newly ghosted entities to the receivers and update the comm list.
  // Unpacking must proceed in entity-rank order so that higher ranking
  // entities that have relations to lower ranking entities will have
  // the lower ranking entities unpacked first.  The higher and lower
  // ranking entities may be owned by different processes,
  // as such unpacking must be performed in rank order.

  //Start of CommAll section:
  {
    stk::CommAll comm( MPI_COMM_WORLD );

    for ( std::set< stk::mesh::EntityProc , stk::mesh::EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {
      stk::mesh::Entity entity = j->first ;
      if ( ! bulk.in_ghost( ghosts , bulk.entity_key(entity) , j->second ) ) {
        // Not already being sent , must send it.
        stk::CommBuffer & buf = comm.send_buffer( j->second );
        buf.pack<unsigned>( bulk.entity_rank(entity) );
        stk::mesh::pack_entity_info(bulk,  buf , entity );
        stk::mesh::pack_field_values(bulk, buf , entity );
      }
    }

    comm.allocate_buffers( size / 4 );

    for ( std::set< stk::mesh::EntityProc , stk::mesh::EntityLess >::iterator
          j = new_send.begin(); j != new_send.end() ; ++j ) {
      stk::mesh::Entity entity = j->first;
      if ( ! bulk.in_ghost( ghosts , bulk.entity_key(entity) , j->second ) ) {
        // Not already being sent , must send it.
        stk::CommBuffer & buf = comm.send_buffer( j->second );
        buf.pack<unsigned>( bulk.entity_rank(entity) );
        stk::mesh::pack_entity_info(bulk,  buf , entity );
        stk::mesh::pack_field_values(bulk, buf , entity );

      }
    }

    comm.communicate();

    std::ostringstream error_msg ;

    for ( int rank = 0 ; rank < rank_count ; ++rank ) {

      for ( int p = 0 ; p < size ; ++p ) {

        stk::CommBuffer & buf = comm.recv_buffer(p);

        while ( buf.remaining() ) {

          // Only unpack if of the current entity rank.
          // If not the current entity rank, break the iteration
          // until a subsequent entity rank iteration.
          {
            int this_rank = ~0u ;
            buf.peek<int>( this_rank );
            if ( this_rank != rank ) break ;

            buf.unpack<int>( this_rank );
          }

          // FIXME for Carol; the code below did not work with -np 4
          //STKUNIT_ASSERT_EQUAL( stk::mesh::unpack_field_values( buf , elem2 , error_msg ), false);
	  //std::cout << "Error message for unpack_field_values = " << error_msg.str() << std::endl ;

        }
      }

    }
  }//end of CommAll section

  bulk.modification_end ();
}

STKUNIT_UNIT_TEST ( UnitTestBulkData_new , testUninitializedMetaData )
{
  stk::ParallelMachine pm = MPI_COMM_WORLD;

  stk::mesh::MetaData meta; // Construct, but do not initialize
  stk::mesh::BulkData bulk(meta, pm);

  meta.initialize(2);

  meta.commit();

  bulk.modification_begin();

  STKUNIT_ASSERT_THROW( bulk.declare_entity(stk::topology::NODE_RANK,
                                            1, /*id*/
                                            stk::mesh::PartVector() ),
                        std::logic_error);
}

namespace {

void new_insert_transitive_closure( stk::mesh::BulkData& bulk_data, std::set<stk::mesh::EntityProc,stk::mesh::EntityLess> & new_send ,
                                const stk::mesh::EntityProc & entry )
{
  // Do not insert if I can determine that this entity is already
  // owned or shared by the receiving processor.

  if ( entry.second != bulk_data.parallel_owner_rank(entry.first) &&
       ! bulk_data.in_shared( bulk_data.entity_key(entry.first), entry.second ) ) {

    std::pair< std::set<stk::mesh::EntityProc,stk::mesh::EntityLess>::iterator , bool >
      result = new_send.insert( entry );

    if ( result.second ) {
      // A new insertion, must also insert the closure

      const unsigned erank = bulk_data.entity_rank(entry.first);

      for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK; irank < erank; ++irank)
      {
        stk::mesh::Entity const *rels_i = bulk_data.begin(entry.first, irank);
        stk::mesh::Entity const *rels_e = bulk_data.end(entry.first, irank);
        for ( ; rels_i != rels_e; ++rels_i)
        {
          stk::mesh::EntityProc tmp( *rels_i , entry.second );
          new_insert_transitive_closure( bulk_data, new_send , tmp );
        }
      }
    }
  }
}


// Synchronize the send list to the receive list.

void new_comm_sync_send_recv(
  stk::mesh::BulkData & mesh ,
  std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > & new_send ,
  std::set< stk::mesh::Entity , stk::mesh::EntityLess > & new_recv )
{
  const int parallel_rank = mesh.parallel_rank();
  const int parallel_size = mesh.parallel_size();

  stk::CommAll all( mesh.parallel() );

  // Communication sizing:

  for ( std::set< stk::mesh::EntityProc , stk::mesh::EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(i->first);
    all.send_buffer( i->second ).skip<stk::mesh::EntityKey>(1).skip<int>(1);
    if ( owner != parallel_rank ) {
      all.send_buffer( owner ).skip<stk::mesh::EntityKey>(1).skip<int>(1);
    }
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  // Communication packing (with message content comments):
  for ( std::set< stk::mesh::EntityProc , stk::mesh::EntityLess >::iterator
        i = new_send.begin() ; i != new_send.end() ; ) {
    const int owner = mesh.parallel_owner_rank(i->first);

    // Inform receiver of ghosting, the receiver does not own
    // and does not share this entity.
    // The ghost either already exists or is a to-be-done new ghost.
    // This status will be resolved on the final communication pass
    // when new ghosts are packed and sent.

    const stk::mesh::EntityKey entity_key = mesh.entity_key(i->first);
    const int proc = i->second;

    all.send_buffer( i->second ).pack(entity_key).pack(proc);

    if ( owner != parallel_rank ) {
      // I am not the owner of this entity.
      // Inform the owner of this ghosting need.
      all.send_buffer( owner ).pack(entity_key).pack(proc);

      // Erase it from my processor's ghosting responsibility:
      // The iterator passed to the erase method will be invalidated.
      std::set< stk::mesh::EntityProc , stk::mesh::EntityLess >::iterator jrem = i ; ++i ;
      new_send.erase( jrem );
    }
    else {
      ++i ;
    }
  }

  all.communicate();

  // Communication unpacking:
  for ( int p = 0 ; p < parallel_size ; ++p ) {
    stk::CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {

      stk::mesh::EntityKey entity_key;
      int proc = 0;

      buf.unpack(entity_key).unpack(proc);

      stk::mesh::Entity const e = mesh.get_entity( entity_key );

      if ( parallel_rank != proc ) {
        //  Receiving a ghosting need for an entity I own.
        //  Add it to my send list.
        STKUNIT_ASSERT( mesh.is_valid(e) );
        stk::mesh::EntityProc tmp( e , proc );
        new_send.insert( tmp );
      }
      else if ( mesh.is_valid(e) ) {
        //  I am the receiver for this ghost.
        //  If I already have it add it to the receive list,
        //  otherwise don't worry about it - I will receive
        //  it in the final new-ghosting communication.
        new_recv.insert( e );
      }
    }
  }
}

void new_comm_recv_to_send(
  stk::mesh::BulkData & mesh ,
  const std::set< stk::mesh::Entity , stk::mesh::EntityLess > & new_recv ,
        std::set< stk::mesh::EntityProc , stk::mesh::EntityLess > & new_send )
{
  const int parallel_size = mesh.parallel_size();

  stk::CommAll all( mesh.parallel() );

  for ( std::set< stk::mesh::Entity , stk::mesh::EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(*i);
    all.send_buffer( owner ).skip<stk::mesh::EntityKey>(1);
  }

  all.allocate_buffers( parallel_size / 4 , false /* Not symmetric */ );

  for ( std::set< stk::mesh::Entity , stk::mesh::EntityLess >::const_iterator
        i = new_recv.begin() ; i != new_recv.end() ; ++i ) {
    const int owner = mesh.parallel_owner_rank(*i);
    const stk::mesh::EntityKey key = mesh.entity_key(*i);
    all.send_buffer( owner ).pack<stk::mesh::EntityKey>( & key , 1 );
  }

  all.communicate();

  for ( int p = 0 ; p < parallel_size ; ++p ) {
    stk::CommBuffer & buf = all.recv_buffer(p);
    while ( buf.remaining() ) {
      stk::mesh::EntityKey key ;
      buf.unpack<stk::mesh::EntityKey>( & key , 1 );
      stk::mesh::EntityProc tmp( mesh.get_entity( key.rank(), key.id() ) , p );
      new_send.insert( tmp );
    }
  }
}

}
