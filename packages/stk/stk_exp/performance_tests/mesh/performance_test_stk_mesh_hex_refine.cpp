#ifndef __IBMCPP__
#include <gtest/gtest.h>
#include <boost/unordered_map.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

#include <performance_tests/mesh/calculate_centroid.hpp>

#include <performance_tests/mesh/hex_refine_info.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <iostream>

#include <boost/range.hpp>

using namespace stk::mesh;

namespace sierra {
namespace mesh {
namespace performance_tests {

namespace {

  void create_entities( BulkData & bulk,
                        Part& node_part,
                        Part& hex_part,
                        // FIXME Part& active_elements_part,
                        HexRefineInfo& refine_info)
  {

    HexRefineInfo refine_info_half(refine_info.m_level-1, refine_info.m_nx, refine_info.m_ny, refine_info.m_nz);
    PartVector add_parts;
    add_parts.push_back(&node_part);

    unsigned eid_start = 1 + refine_info.elem_id_offset();
    unsigned eid_end = eid_start + refine_info.num_elems();

    unsigned nid_start = 1 + refine_info.node_id_offset();
    unsigned nid_end = nid_start + refine_info.num_nodes();

    std::cout << "nid_start = " << nid_start << " nid_end= " << nid_end << " diff= " << nid_end - nid_start << std::endl;
    std::cout << "eid_start = " << eid_start << " eid_end= " << eid_end << " diff= " << eid_end - eid_start << std::endl;

    boost::unordered_map<unsigned, Entity*> node_map;

    for(unsigned nid=nid_start; nid<nid_end; ++nid) {
      Entity *key = &bulk.declare_entity(0, nid, add_parts);
      node_map[nid] = key;
    }

    for (unsigned entity_id = eid_start; entity_id < eid_end; ++entity_id)  {
      unsigned ix = 0, iy = 0, iz = 0;
      refine_info.elem_x_y_z(entity_id, ix, iy, iz);
      EntityId ie_check = refine_info.elem_id(ix, iy, iz);
      EXPECT_EQ(ie_check, entity_id);

      stk::mesh::EntityId elem_node[8] ;

      elem_node[0] = refine_info.node_id( ix   , iy   , iz   );
      elem_node[1] = refine_info.node_id( ix+1 , iy   , iz   );
      elem_node[2] = refine_info.node_id( ix+1 , iy   , iz+1 );
      elem_node[3] = refine_info.node_id( ix   , iy   , iz+1 );
      elem_node[4] = refine_info.node_id( ix   , iy+1 , iz   );
      elem_node[5] = refine_info.node_id( ix+1 , iy+1 , iz   );
      elem_node[6] = refine_info.node_id( ix+1 , iy+1 , iz+1 );
      elem_node[7] = refine_info.node_id( ix   , iy+1 , iz+1 );

      // check if a parent node
      for (unsigned i = 0; i<8; ++i) {
        unsigned ixn = 0, iyn = 0, izn = 0;
        refine_info.node_x_y_z(elem_node[i], ixn, iyn, izn);
        EntityId in_check = refine_info.node_id(ixn, iyn, izn);
        EXPECT_EQ(in_check, elem_node[i]);

        if (
            ((ixn - 1) % 2 == 0) &&
            ((iyn - 1) % 2 == 0) &&
            ((izn - 1) % 2 == 0))
          {
            elem_node[i] = refine_info_half.node_id(ixn/2, iyn/2, izn/2);
          }
      }

      stk::mesh::declare_element( bulk, hex_part, entity_id , elem_node);

#if 0
      for (unsigned i = 0; i<8; ++i) {
        stk::mesh::Entity * const node = bulk.get_entity( MetaData::NODE_RANK , elem_node[i] );
        bulk.change_entity_parts(*node, add_parts);

        ThrowRequireMsg( node != NULL, "found null node in create_entities");

        // Compute and assign coordinates to the node
        unsigned nx = 0, ny = 0, nz = 0;
        refine_info.node_x_y_z(elem_node[i], nx, ny, nz);

        Scalar * data = stk::mesh::field_data( m_coord_field , *node );

        data[0] = nx ;
        data[1] = ny ;
        data[2] = -(Scalar)nz ;
      }
#endif
    }
  }
}

TEST( Hex, STKMesh_Performance_on_Hex_Refine)
{
  double start_time = stk::cpu_time();
  //unsigned ex=100, ey=100, ez=100;
  const int max_levels = 2;
  unsigned nn = 50/(1 << (max_levels-1));
  unsigned ex=nn, ey=nn, ez=nn;

  //unsigned num_elems = ex*ey*ez;
  fixtures::HexFixture fixture(MPI_COMM_SELF, ex, ey, ez);
  fixture.m_fem_meta.commit();
  fixture.generate_mesh();
  double end_time = stk::cpu_time() - start_time;

  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3, 0.0);
  //const double tolerance = 1.e-6;
  //const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  for (int level=1; level <= max_levels; ++level) {

    HexRefineInfo refine_info(level, ex, ey, ez);

    //Selector hex_elem_selector(fixture.m_hex_part & !fixture.m_node_part);

    unsigned num_elems_new = refine_info.num_elems();
    std::cout << "num_elems_new for level = " << level << " = " << num_elems_new << std::endl;
    fixture.m_bulk_data.modification_begin();
    create_entities(fixture.m_bulk_data, fixture.m_node_part, fixture.m_hex_part, refine_info);
    fixture.m_bulk_data.modification_end();
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Num refines: " << max_levels << std::endl;
  std::cout << "Time to refine: " << test_time << std::endl;
}

}
}
}
#endif
