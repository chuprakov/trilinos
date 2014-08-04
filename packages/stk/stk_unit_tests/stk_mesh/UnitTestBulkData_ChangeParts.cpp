/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>                     // for size_t
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <vector>                       // for vector
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Types.hpp"      // for PartVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Selector; } }

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Entity;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::fixtures::RingFixture;

//----------------------------------------------------------------------

TEST(UnitTestingOfBulkData, testChangeParts_ringmesh)
{
  // This unit test tests part operations and verifies operations
  // by looking at bucket supersets.

  stk::ParallelMachine pm = MPI_COMM_WORLD;
  MPI_Barrier( pm );

  const unsigned nPerProc   = 10;
  const int p_rank     = stk::parallel_machine_rank( pm );
  const int p_size     = stk::parallel_machine_size( pm );
  const unsigned nLocalNode = nPerProc + ( 1 < p_size ? 1 : 0 );
  const unsigned nLocalElement = nPerProc ;

  // Create the ring mesh

  RingFixture ring_mesh( pm , nPerProc , true /* generate parts */ );
  ring_mesh.m_meta_data.commit();
  BulkData& bulk = ring_mesh.m_bulk_data;

  bulk.modification_begin();
  ring_mesh.generate_mesh( );
  ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, false /* no aura */));

  bulk.modification_begin();
  ring_mesh.fixup_node_ownership();
  ASSERT_TRUE(stk::unit_test::modification_end_wrapper(bulk, false /* no aura */));

  Part & part_owns = ring_mesh.m_meta_data.locally_owned_part();
  Part & part_univ = ring_mesh.m_meta_data.universal_part();

  // Check that local elements are in the expected parts. Note that the
  // RingMesh puts each element in its own part.
  for ( unsigned i = 0 ; i < nLocalElement ; ++i ) {
    const unsigned n = i + nPerProc * p_rank ;
    Entity const element = bulk.get_entity( stk::topology::ELEMENT_RANK /*entity rank*/,
                                              ring_mesh.m_element_ids[n] );
    ASSERT_TRUE( bulk.is_valid(element) );
    ASSERT_TRUE( bulk.bucket(element).member( part_univ ) );
    ASSERT_TRUE( bulk.bucket(element).member( part_owns ) );
    ASSERT_TRUE( bulk.bucket(element).member( * ring_mesh.m_element_parts[ n % ring_mesh.m_element_parts.size() ] ) );
  }

  // Check that local nodes are in the expected parts. Note that the relations
  // that nodes have to elements should cause induced membership of the node
  // in the parts of both elements it touches.
  for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
    const unsigned n = ( i + nPerProc * p_rank ) % ring_mesh.m_node_ids.size();
    const unsigned e0 = n ;
    const unsigned e1 = ( n + ring_mesh.m_element_ids.size() - 1 ) % ring_mesh.m_element_ids.size();
    const unsigned ns = ring_mesh.m_element_parts.size();
    const unsigned n0 = e0 % ns ;
    const unsigned n1 = e1 % ns ;
    Part * const epart_0 = ring_mesh.m_element_parts[ n0 < n1 ? n0 : n1 ];
    Part * const epart_1 = ring_mesh.m_element_parts[ n0 < n1 ? n1 : n0 ];

    Entity const node = bulk.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[n] );
    ASSERT_TRUE( bulk.is_valid(node) );
    if ( bulk.parallel_owner_rank(node) == p_rank ) {
      ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
      ASSERT_TRUE( bulk.bucket(node).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(node).member( *epart_0 ) );
      ASSERT_TRUE( bulk.bucket(node).member( *epart_1 ) );
    }
    else {
      ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
      ASSERT_TRUE( ! bulk.bucket(node).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(node).member( * epart_0 ) );
      ASSERT_TRUE( bulk.bucket(node).member( * epart_1 ) );
    }
  }

  bulk.modification_begin();

  // On rank 0, change all locally owned elements to the extra-part then check
  // for correct part membership
  if ( 0 == p_rank ) {
    for ( unsigned i = 0 ; i < nLocalElement ; ++i ) {
      const unsigned n = i + nPerProc * p_rank ;

      PartVector add(1); add[0] = & ring_mesh.m_element_part_extra ;
      PartVector rem(1); rem[0] = ring_mesh.m_element_parts[ n % ring_mesh.m_element_parts.size() ];

      Entity const element = bulk.get_entity( stk::topology::ELEMENT_RANK , ring_mesh.m_element_ids[n] );
      bulk.change_entity_parts( element , add , rem );
      ASSERT_TRUE( bulk.bucket(element).member( part_univ ) );
      ASSERT_TRUE( bulk.bucket(element).member( part_owns ) );
      ASSERT_TRUE( bulk.bucket(element).member(ring_mesh.m_element_part_extra ) );
    }
  }

  bulk.modification_end();

  // Modification end has been called, check that the part changes made
  // in the previous step are reflected across the other procs.
  for ( unsigned i = 0 ; i < nLocalNode ; ++i ) {
    const unsigned n = ( i + nPerProc * p_rank ) % ring_mesh.m_node_ids.size();
    const unsigned e0 = n ;
    const unsigned e1 = ( n + ring_mesh.m_element_ids.size() - 1 ) % ring_mesh.m_element_ids.size();
    const unsigned ns = ring_mesh.m_element_parts.size();
    const unsigned n0 = e0 % ns ;
    const unsigned n1 = e1 % ns ;
    Part * ep_0 = e0 < nLocalElement ? & ring_mesh.m_element_part_extra : ring_mesh.m_element_parts[n0] ;
    Part * ep_1 = e1 < nLocalElement ? & ring_mesh.m_element_part_extra : ring_mesh.m_element_parts[n1] ;

    Part * epart_0 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_0 : ep_1 ;
    Part * epart_1 = ep_0->mesh_meta_data_ordinal() < ep_1->mesh_meta_data_ordinal() ? ep_1 : ep_0 ;

    Entity const node = bulk.get_entity( stk::topology::NODE_RANK , ring_mesh.m_node_ids[n] );
    ASSERT_TRUE( bulk.is_valid(node) );
    if ( bulk.parallel_owner_rank(node) == p_rank ) {
      ASSERT_TRUE( bulk.bucket(node).member( part_owns ) );
    }
    else {
      ASSERT_TRUE( ! bulk.bucket(node).member( part_owns ) );
    }

#if 0
    //// DEBUG

    std::cout << "p_rank = " << p_rank << "; i = " << i << std::endl;
    std::cout << " node of interest " << node << std::endl;
#ifdef USE_STK_MESH_IMPL_PARTITION
    std::cout << " its partition " << *node.bucket().getPartition() << std::endl;
#else
    std::cout << " its bucket " << node.bucket() << std::endl;
#endif
    std::cout << " epart_0 " << *epart_0 << " " << epart_0->mesh_meta_data_ordinal() << std::endl;
    std::cout << " epart_1 " << *epart_1 << " " << epart_1->mesh_meta_data_ordinal() << std::endl;

    //// GUBED
#endif

    ASSERT_TRUE( bulk.bucket(node).member( part_univ ) );
    ASSERT_TRUE( bulk.bucket(node).member( *epart_0 ) );
    ASSERT_TRUE( bulk.bucket(node).member( *epart_1 ) );
  }
}
