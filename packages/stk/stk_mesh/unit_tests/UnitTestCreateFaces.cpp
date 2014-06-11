/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
/*
 * UnitTestCreateFaces.C created by tcfishe on Feb 20, 2014
 */

#include <stddef.h>                     // for size_t
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Comm.hpp>       // for comm_mesh_counts
#include <stk_mesh/base/CreateFaces.hpp> // for create_faces
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/SkinMesh.hpp>   // for skin_mesh
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture
#include <stk_mesh/fixtures/TetFixture.hpp>  // for TetFixture
#include <stk_mesh/fixtures/GearsFixture.hpp>  // for GearsFixture
#include <stk_mesh/fixtures/heterogeneous_mesh.hpp>
#include <stk_mesh/fixtures/degenerate_mesh.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <vector>                       // for vector, vector<>::iterator
#include "gtest/gtest.h"                // for AssertHelper
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/Types.hpp"      // for BucketVector, EntityRank
#include "stk_topology/topology.hpp"    // for topology, etc

using stk::mesh::MetaData;

namespace {
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
  const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

  const size_t nodes_per_hex = 8;
  const size_t faces_per_hex = 6;
  const size_t nodes_per_quad= 4;
  
  const size_t nodes_per_tet = 4;
  const size_t faces_per_tet = 4;
  const size_t nodes_per_tri = 3;
  
  size_t exp_hex_face_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_face = 3 * nx * ny * nz;
    exp_face += ny * nz + nz * nx + nx * ny;
    return exp_face;
  }

  size_t exp_tet_face_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_face = 12 * nx * ny * nz;
    exp_face += 2* (ny * nz + nz * nx + nx * ny);
    return exp_face;
  }

  size_t exp_node_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_node = (nx+1) * (ny+1) * (nz+1);
    return exp_node;
  }

  size_t exp_hex_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_elem = nx * ny * nz;
    return exp_elem;
  }

  size_t exp_tet_count(size_t nx, size_t ny, size_t nz)
  {
    size_t exp_elem = 6 * nx * ny * nz;
    return exp_elem;
  }
}

STKUNIT_UNIT_TEST ( UnitTestCreateFaces, Hex_2x1x1 )
{
  const size_t NX = 2;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

}

STKUNIT_UNIT_TEST ( UnitTestCreateFaces, Tet_2x1x1 )
{
  const size_t NX = 2;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

}

STKUNIT_UNIT_TEST( UnitTestCreateFaces , Hex_3x1x1 )
{
  const size_t NX = 3;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateFaces , Tet_3x1x1 )
{
  const size_t NX = 3;
  const size_t NY = 1;
  const size_t NZ = 1;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }
}

STKUNIT_UNIT_TEST( UnitTestCreateFaces , testCreateFaces3x3x3 )
{
  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( faces_per_hex, b.num_faces(i) );
      STKUNIT_EXPECT_EQ( 0u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( nodes_per_hex,  b.num_nodes(i) );
    }
  }

  stk::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = face_buckets.begin();
       b_itr != face_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( 0u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( nodes_per_quad, b.num_nodes(i) );
    }
  }

  unsigned num_ghosted_faces = 0;

  stk::mesh::Selector ghosted_part = fixture.m_meta.universal_part() &
      !(fixture.m_meta.locally_owned_part() | fixture.m_meta.globally_shared_part());

  stk::mesh::BucketVector  ghosted_elem_buckets = fixture.m_bulk_data.get_buckets(elem_rank,ghosted_part);

  for ( stk::mesh::BucketVector::iterator b_itr = ghosted_elem_buckets.begin();
      b_itr != ghosted_elem_buckets.end();
      ++b_itr
  )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      const stk::mesh::Entity & elem = b[i];

      const unsigned num_faces = fixture.m_bulk_data.num_faces(elem);
      STKUNIT_EXPECT_EQ(faces_per_hex, num_faces);

      const stk::mesh::Entity * faces = fixture.m_bulk_data.begin_faces(elem);
      unsigned elem_num_ghosted_faces = 0;
      unsigned elem_num_shared_faces = 0;

      const int my_proc_id = fixture.m_bulk_data.parallel_rank();

      for (unsigned k = 0; k < num_faces; ++k)
      {
        const stk::mesh::Entity & face = faces[k];

        const bool is_ghosted = fixture.m_bulk_data.in_receive_ghost(fixture.m_bulk_data.entity_key(face));
        if (is_ghosted)
        {
          elem_num_ghosted_faces += 1;
          const int face_owner = fixture.m_bulk_data.parallel_owner_rank(face);
          STKUNIT_EXPECT_NE(face_owner, my_proc_id);
        }

        const bool is_shared = fixture.m_bulk_data.in_shared(fixture.m_bulk_data.entity_key(face));
        if (is_shared) elem_num_shared_faces += 1;

        const bool is_ghosted_or_shared = is_ghosted || is_shared;
        STKUNIT_ASSERT(is_ghosted_or_shared);

      }

      STKUNIT_EXPECT_LE(1u, elem_num_ghosted_faces);
      num_ghosted_faces += elem_num_ghosted_faces;

      const unsigned shared_and_ghosted = elem_num_shared_faces + elem_num_ghosted_faces;
      STKUNIT_EXPECT_EQ(faces_per_hex, shared_and_ghosted);
    }
  }

  if (fixture.m_bulk_data.parallel_size() > 1)
    STKUNIT_EXPECT_LE(1u, num_ghosted_faces);


}

STKUNIT_UNIT_TEST( UnitTestCreateFaces , testCreateTetFaces3x3x3 )
{
  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::TetFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_tet_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_tet_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( faces_per_tet, b.num_faces(i) );
      STKUNIT_EXPECT_EQ( 0u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( nodes_per_tet,  b.num_nodes(i) );
    }
  }

  stk::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = face_buckets.begin();
       b_itr != face_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      STKUNIT_EXPECT_EQ( 0u, b.num_edges(i) );
      STKUNIT_EXPECT_EQ( nodes_per_tri, b.num_nodes(i) );
    }
  }

  unsigned num_ghosted_faces = 0;

  stk::mesh::Selector ghosted_part = fixture.m_meta.universal_part() &
      !(fixture.m_meta.locally_owned_part() | fixture.m_meta.globally_shared_part());

  stk::mesh::BucketVector  ghosted_elem_buckets = fixture.m_bulk_data.get_buckets(elem_rank,ghosted_part);

  for ( stk::mesh::BucketVector::iterator b_itr = ghosted_elem_buckets.begin();
      b_itr != ghosted_elem_buckets.end();
      ++b_itr
  )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      const stk::mesh::Entity & elem = b[i];

      const unsigned num_faces = fixture.m_bulk_data.num_faces(elem);
      STKUNIT_EXPECT_EQ(faces_per_tet, num_faces);

      const stk::mesh::Entity * faces = fixture.m_bulk_data.begin_faces(elem);
      unsigned elem_num_ghosted_faces = 0;
      unsigned elem_num_shared_faces = 0;

      const int my_proc_id = fixture.m_bulk_data.parallel_rank();

      for (unsigned k = 0; k < num_faces; ++k)
      {
        const stk::mesh::Entity & face = faces[k];

        const bool is_ghosted = fixture.m_bulk_data.in_receive_ghost(fixture.m_bulk_data.entity_key(face));
        if (is_ghosted)
        {
          elem_num_ghosted_faces += 1;
          const int face_owner = fixture.m_bulk_data.parallel_owner_rank(face);
          STKUNIT_EXPECT_NE(face_owner, my_proc_id);
        }

        const bool is_shared = fixture.m_bulk_data.in_shared(fixture.m_bulk_data.entity_key(face));
        if (is_shared) elem_num_shared_faces += 1;

        const bool is_ghosted_or_shared = is_ghosted || is_shared;
        STKUNIT_ASSERT(is_ghosted_or_shared);

      }

      STKUNIT_EXPECT_LE(1u, elem_num_ghosted_faces);
      num_ghosted_faces += elem_num_ghosted_faces;

      const unsigned shared_and_ghosted = elem_num_shared_faces + elem_num_ghosted_faces;
      STKUNIT_EXPECT_EQ(faces_per_tet, shared_and_ghosted);
    }
  }

  if (fixture.m_bulk_data.parallel_size() > 1)
    STKUNIT_EXPECT_LE(1u, num_ghosted_faces);


}

STKUNIT_UNIT_TEST( UnitTestCreateFaces , testSkinAndCreateFaces3x3x3 )
{
  const stk::mesh::EntityRank elem_rank = stk::topology::ELEMENT_RANK;
  const stk::mesh::EntityRank face_rank = stk::topology::FACE_RANK;
  const stk::mesh::EntityRank edge_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityRank node_rank = stk::topology::NODE_RANK;

  const size_t NX = 3;
  const size_t NY = 3;
  const size_t NZ = 3;

  stk::mesh::fixtures::HexFixture fixture( MPI_COMM_WORLD, NX, NY, NZ);

  fixture.m_meta.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.m_bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.m_bulk_data , counts);

    STKUNIT_EXPECT_EQ( exp_node_count(NX, NY, NZ), counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( exp_hex_face_count(NX, NY, NZ), counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( exp_hex_count(NX, NY, NZ), counts[elem_rank] ); // elements
  }

  stk::mesh::BucketVector  elem_buckets = fixture.m_bulk_data.buckets(elem_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = elem_buckets.begin();
       b_itr != elem_buckets.end();
       ++b_itr
      )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned elem_ordinal = i;
      STKUNIT_EXPECT_EQ( faces_per_hex,  b.num_faces(elem_ordinal) );
      STKUNIT_EXPECT_EQ( nodes_per_hex,  b.num_nodes(elem_ordinal) );
    }
  }

  stk::mesh::BucketVector  face_buckets = fixture.m_bulk_data.buckets(face_rank);
  for ( stk::mesh::BucketVector::iterator b_itr = face_buckets.begin();
      b_itr != face_buckets.end();
      ++b_itr
  )
  {
    stk::mesh::Bucket & b = **b_itr;
    for ( size_t i = 0; i< b.size(); ++i) {
      unsigned face_ordinal = i;
      STKUNIT_EXPECT_EQ( 0u, b.num_edges(face_ordinal) );
      STKUNIT_EXPECT_EQ( nodes_per_quad, b.num_nodes(face_ordinal) );
    }
  }

}

STKUNIT_UNIT_TEST ( UnitTestCreateFaces, Gears )
{
  stk::mesh::fixtures::GearsFixture fixture( MPI_COMM_WORLD, 1,
					     stk::mesh::fixtures::GearParams(0.1, 0.4, 1.0, -0.4, 0.4));

  fixture.meta_data.commit();
  fixture.generate_mesh();

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.bulk_data , counts);

    STKUNIT_EXPECT_EQ( 4340u,  counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,                         counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,                         counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 3348u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(fixture.bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts( fixture.bulk_data , counts);

    STKUNIT_EXPECT_EQ( 4340u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u                        , counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 10974u, counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 3348u, counts[elem_rank] ); // elements
  }

}

STKUNIT_UNIT_TEST ( UnitTestCreateFaces, Heterogeneous )
{
  stk::mesh::MetaData meta_data(3);
  stk::mesh::fixtures::VectorFieldType & node_coord =
    meta_data.declare_field<stk::mesh::fixtures::VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field( node_coord , meta_data.universal_part() , 3);

  stk::mesh::fixtures::heterogeneous_mesh_meta_data( meta_data , node_coord );
  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, MPI_COMM_WORLD );
  stk::mesh::fixtures::heterogeneous_mesh_bulk_data( bulk_data , node_coord );

  /*
   * Three hexes, three wedges, three tets, two pyramids,
   * three quad shells, and three triangle shells.
   *
   *  Z = 0 plane:
   *
   *    Y
   *    ^   9      10
   *    !   *-------*
   *    !  / \     / \
   *    ! /   \   /   \
   *     /     \ /     \
   *    *-------*-------*-------*
   *   5|      6|      7|      8|
   *    |       |       |       |
   *    |       |       |       |
   *    *-------*-------*-------*    ----> X
   *    1       2       3       4
   *
   *  Z = -1 plane:
   *
   *    Y
   *    ^  19      20
   *    !   *-------*
   *    !  / \     / \
   *    ! /   \   /   \
   *     /     \ /     \
   *    *-------*-------*-------*
   *  15|     16|     17|     18|
   *    |       |       |       |
   *    |       |       |       |
   *    *-------*-------*-------*    ----> X
   *   11      12      13      14
   *
   *
   *  Last node (#21) at Z = -2, translated from node #16
   */

  /*
   * Total elements = 21 = 3 + 3 + 3 + 2 + 3 + 3
   * Total faces    = 39 = 6 front + 6 back + 9 perimeter + 6 internal + 7 from tet + 5 from pyr
   */
  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    STKUNIT_EXPECT_EQ( 21u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,  counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u,  counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 17u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    STKUNIT_EXPECT_EQ( 21u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u,  counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 39u, counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 17u, counts[elem_rank] ); // elements
  }

}

STKUNIT_UNIT_TEST ( UnitTestCreateFaces, Degenerate )
{
  stk::mesh::MetaData meta_data(3);
  stk::mesh::fixtures::VectorFieldType & node_coord =
    meta_data.declare_field<stk::mesh::fixtures::VectorFieldType>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::put_field( node_coord , meta_data.universal_part() , 3);

  stk::mesh::fixtures::degenerate_mesh_meta_data( meta_data , node_coord );
  meta_data.commit();

  stk::mesh::BulkData bulk_data( meta_data, MPI_COMM_WORLD );
  stk::mesh::fixtures::degenerate_mesh_bulk_data( bulk_data , node_coord );

  /*
   *  Z = 0 plane:
   *
   *    Y
   *    ^ 
   *    !
   *     
   *   4*       *5
   *    |\     /|
   *    | \   / |
   *    |  \ /  |
   *    *---*---* ----> X
   *    1   2   3
   *
   *  Z = -1 plane:
   *
   *    Y
   *    ^ 
   *    !
   *     
   *   9*       *10
   *    |\     /|
   *    | \   / |
   *    |  \ /  |
   *    *---*---* ----> X
   *    6   7   8
   *
   * Total elements = 2 hexes
   * Total faces    = 11 = 2 front + 2 back + 6 perimeter + 1 internal degenerate
   */
  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    STKUNIT_EXPECT_EQ(10u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u, counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ( 0u, counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 2u, counts[elem_rank] ); // elements
  }

  stk::mesh::create_faces(bulk_data);

  {
    std::vector<size_t> counts ;
    stk::mesh::comm_mesh_counts(bulk_data , counts);

    STKUNIT_EXPECT_EQ(10u, counts[node_rank] ); // nodes
    STKUNIT_EXPECT_EQ( 0u, counts[edge_rank] ); // edges
    STKUNIT_EXPECT_EQ(11u, counts[face_rank] ); // faces
    STKUNIT_EXPECT_EQ( 2u, counts[elem_rank] ); // elements
  }

}

