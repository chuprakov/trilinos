/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for size_t
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <gtest/gtest.h>
#include <string>                       // for string, char_traits
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityId, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { namespace fixtures { class BoxFixture; } } }
namespace stk { namespace mesh { struct EntityKey; } }


using stk::mesh::Part;
using stk::mesh::Bucket;
using stk::mesh::PairIterRelation;
using stk::mesh::PairIterEntityComm;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::BaseEntityRank;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::BucketVector;
using stk::mesh::fixtures::BoxFixture;

namespace {
const EntityRank NODE_RANK = stk::topology::NODE_RANK;
} // empty namespace


void printBuckets(std::ostringstream& msg, BulkData& mesh)
{
  const BucketVector & buckets = mesh.buckets(NODE_RANK);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      const Bucket& bucket = *buckets[i];
      msg << " bucket[" << i << "] = ";
      size_t bucket_size = bucket.size();
      for (unsigned ie=0; ie < bucket_size; ie++)
        {
          msg << mesh.identifier(bucket[ie]) << ", ";
        }
    }
}

static void checkBuckets( BulkData& mesh)
{
  const BucketVector & buckets = mesh.buckets(NODE_RANK);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      Bucket* bucket = buckets[i];
      ASSERT_TRUE(bucket->assert_correct());
    }
}

TEST(UnitTestingOfBulkData, test_other_ghosting_2)
{
  //
  // testing if modification flags propagate properly for ghosted entities
  //
  // To test this, we focus on a single node shared on 2 procs, ghosted on others
  //

  /**
   * 1D Mesh (node,owner)--[elem,owner]---(...)
   *
   * <---(70,0)--[500,0]--(41,1)--[301,1]---(42,2)---[402,2]---(70,0)--->
   *
   * <---(50,0)--[100,0]--(21,1)--[201,1]---(32,2)---[302,2]---(50,0)--->
   *
   */

  // elem, node0, node1, owner
  EntityId elems_0[][4] = { {100, 21, 50, 0}, {201, 21, 32, 1}, {302, 32, 50, 2},
                            {500, 41, 70, 0}, {301, 41, 42, 1}, {402, 42, 70, 2}  };
  // node, owner
  EntityId nodes_0[][2] = { {21,1}, {50,0}, {32, 2}, {41, 1}, {42, 1}, {70, 0} };

  unsigned nelems = sizeof(elems_0)/4/sizeof(EntityId);
  unsigned nnodes = sizeof(nodes_0)/2/sizeof(EntityId);

  stk::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;

  std::vector<std::string> entity_rank_names = stk::mesh::entity_rank_names();
  entity_rank_names.push_back("FAMILY_TREE");

  MetaData meta_data(spatial_dim, entity_rank_names);
  //Part & part_tmp = meta_data.declare_part( "temp");

  meta_data.commit();
  unsigned max_bucket_size = 1;
  BulkData mesh(meta_data, pm, max_bucket_size);
  //BulkData mesh(MetaData::get_meta_data(meta_data), pm);
  int p_rank = mesh.parallel_rank();
  int p_size = mesh.parallel_size();

  if (p_size != 3) return;

  //
  // Begin modification cycle so we can create the entities and relations
  //

  // We're just going to add everything to the universal part
  stk::mesh::PartVector empty_parts;

  // Create elements
  const EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  Entity elem = Entity();

  mesh.modification_begin();

  for (unsigned ielem=0; ielem < nelems; ielem++)
    {
      if (static_cast<int>(elems_0[ielem][3]) == p_rank)
        {
          elem = mesh.declare_entity(elem_rank, elems_0[ielem][0], empty_parts);

          EntityVector nodes;
          // Create node on all procs
          nodes.push_back( mesh.declare_entity(NODE_RANK, elems_0[ielem][2], empty_parts) );
          nodes.push_back( mesh.declare_entity(NODE_RANK, elems_0[ielem][1], empty_parts) );

          // Add relations to nodes
          mesh.declare_relation( elem, nodes[0], 0 );
          mesh.declare_relation( elem, nodes[1], 1 );

        }
    }

  mesh.modification_end();

  Entity node1 = Entity();

  // change node owners
  mesh.modification_begin();

  std::vector<EntityProc> change;

  for (unsigned inode=0; inode < nnodes; inode++)
    {
      node1 = mesh.get_entity(stk::topology::NODE_RANK, nodes_0[inode][0]);
      if (mesh.is_valid(node1) && mesh.parallel_owner_rank(node1) == p_rank)
        {
          int dest = nodes_0[inode][1];
          EntityProc eproc(node1, dest);
          change.push_back(eproc);
        }
    }

  mesh.change_entity_owner( change );

  mesh.modification_end();

  checkBuckets(mesh);

  MPI_Barrier(MPI_COMM_WORLD);


  // attempt to delete a node and its elems but on a ghosted proc
  mesh.modification_begin();

  if (p_rank == 2)
    {
      node1 = mesh.get_entity(stk::topology::NODE_RANK, 21);
      Entity elem1 = mesh.get_entity(stk::topology::ELEMENT_RANK, 201);
      Entity elem2 = mesh.get_entity(stk::topology::ELEMENT_RANK, 100);

      bool did_it_elem = mesh.destroy_entity(elem1);
      did_it_elem = did_it_elem & mesh.destroy_entity(elem2);
      ASSERT_TRUE(did_it_elem);
      bool did_it = mesh.destroy_entity(node1);
      ASSERT_TRUE(did_it);
    }

  mesh.modification_end();

  checkBuckets(mesh);

  // this node should no longer exist anywhere
  node1 = mesh.get_entity(stk::topology::NODE_RANK, 21);

  // uncomment to force failure of test
  // ASSERT_TRUE(node1 == 0);

}

