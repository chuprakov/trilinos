#include <gtest/gtest.h>

#include <sierra/io/array_mesh_fixture.hpp>
#include <sierra/mesh/fixture/array_mesh_hex_fixture.hpp>

#include <stk_util/environment/CPUTime.hpp>

#include <performance_tests/mesh/calculate_centroid.hpp>

#include <cmath>
#include <iostream>

using namespace sierra::mesh;

namespace {

const unsigned spatial_dim = 3;

template<typename Topology>
void block_gather(int num_elems,
            const std::vector<int>& connectivity,
            const std::vector<double>& coord_field,
            std::vector<double>& avg_centroid)
{
  double elem_centroid[spatial_dim] = {0, 0, 0};
  const int num_nodes_per_elem = sierra::mesh::num_nodes<Topology>::value;
  double elem_node_coords[num_nodes_per_elem*spatial_dim];

  //gather nodal coordinates for each element, call calculate_centroid.
  int elem_offset = 0;
  for( int i=0; i<num_elems; ++i) {
    int offset = 0;
    for( int n=0; n<num_nodes_per_elem; ++n) {
      const double * node_coords = &coord_field[connectivity[elem_offset]*spatial_dim];
      elem_offset++;
      elem_node_coords[offset++] = node_coords[0];
      elem_node_coords[offset++] = node_coords[1];
      elem_node_coords[offset++] = node_coords[2];
    }

    performance_tests::calculate_centroid_3d(num_nodes_per_elem,&elem_node_coords[0], elem_centroid);

    //add this element-centroid to the avg_centroid vector, and
    //re-zero the element-centroid vector:
    avg_centroid[0] += elem_centroid[0]; elem_centroid[0] = 0;
    avg_centroid[1] += elem_centroid[1]; elem_centroid[1] = 0;
    avg_centroid[2] += elem_centroid[2]; elem_centroid[2] = 0;
  }
}

void array_mesh_gather( const array_mesh & mesh,
                  const std::vector<double>& coord_field,
                  std::vector<double> & avg_centroid )
{
  array_mesh::BlockRange blocks = mesh.get_blocks();

  for(array_mesh::BlockIterator b_it=blocks.first, b_end=blocks.second; b_it!=b_end; ++b_it) {
    int num_elems = mesh.get_num_elems(*b_it);
    int topology = mesh.get_topology(*b_it);

    switch(topology) {
     case Hex8::value:
       block_gather<Hex8>(num_elems, mesh.get_block_connectivity(*b_it),
                            coord_field, avg_centroid);
       break;
     case Node::value: continue; break;
     case Tet4::value:
     default:
       std::cout<<"Unsupported topology!"<<std::endl;
    }
  }
}

}

TEST( ArrayMesh, gather_centroid_hex_elem_genmesh)
{
  double start_time = stk::cpu_time();
#ifndef NDEBUG
  unsigned ex=2,  ey=2,  ez=2; // make things smaller in debug
#else
  unsigned ex=100, ey=100, ez=100;
#endif
  unsigned num_elems = ex*ey*ez;

  std::ostringstream oss;
  oss << ex<<"x"<<ey<<"x"<<ez;
  std::string file_name = oss.str();
  std::string mesh_type("generated");

  sierra::mesh::io::array_mesh_fixture fixture(MPI_COMM_WORLD, mesh_type, file_name);

  double end_time = stk::cpu_time() - start_time;

  std::cout << "\tNum Nodes: " << fixture.m_mesh.get_num_nodes()<<std::endl;
  std::cout << "\tNum Elements: " << fixture.m_mesh.get_num_elements() << std::endl;
  std::cout << "Time to create hex mesh: " << end_time << std::endl;

  std::vector<double> avg_centroid(3,0);
  const double tolerance = 1.e-6;
  const double expected = ((double)ex)/2;

  start_time = stk::cpu_time();

  const int spatial_dim = 3;
  const int num_iters = 100;
  for(int t=0; t<num_iters; ++t) {
    array_mesh_gather(fixture.m_mesh, fixture.m_coords, avg_centroid);

    for(int d=0; d<spatial_dim; ++d) {
      EXPECT_LT(std::abs(avg_centroid[d]/num_elems - expected), tolerance);
    }

    avg_centroid[0] = 0;
    avg_centroid[1] = 0;
    avg_centroid[2] = 0;
  }

  end_time = stk::cpu_time();

  double test_time = end_time - start_time;
  std::cout << "Time to compute centroids ("<<num_iters<<" iters): " << test_time << std::endl;
}
