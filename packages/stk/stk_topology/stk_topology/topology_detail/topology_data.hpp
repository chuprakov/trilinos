#ifndef STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP
#define STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP

#include <stk_topology/topology.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/vector/vector30_c.hpp>
#include <boost/utility.hpp>

//TODO implement permutations for tets, pyramids, wedges and hexes
//TODO implement permutations polarity

namespace stk { namespace topology_detail {

template <topology::topology_t Topology, typename Enable = void>
struct topology_data;


//***************************************************************************
// topology::INVALID -- topology::INVALID_TOPOLOGY
//***************************************************************************

template <>
struct topology_data<topology::INVALID_TOPOLOGY>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::INVALID_TOPOLOGY;
  static const topology::topology_t base = topology::INVALID_TOPOLOGY;

  static const bool is_valid = false;
  static const topology::rank_t rank = topology::INVALID_RANK;
  static const topology::rank_t side_rank = topology::INVALID_RANK;
  static const topology::topology_t edge_topology = topology::INVALID_TOPOLOGY;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 0;
  static const unsigned num_nodes = 0;
  static const unsigned num_vertices = 0;
  static const unsigned num_edges = 0;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 0;
  static const unsigned num_positive_permutations = 0;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , false // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<> edge_node_ordinals_vector;
  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<> permutation_node_ordinals_vector;
};

//***************************************************************************
// topology::HETEROGENEOUS_EDGE -- topology::HETEROGENEOUS_EDGE
//***************************************************************************

template <>
struct topology_data<topology::HETEROGENEOUS_EDGE>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static const topology::topology_t value = topology::HETEROGENEOUS_EDGE;
  static const topology::topology_t base = value;
  static const bool is_valid = true;
  static const topology::rank_t rank = topology::EDGE_RANK;
  static const topology::rank_t side_rank = topology::NODE_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true // 2d
                                , true // 3d
                              > spatial_dimension_vector;

};

//***************************************************************************
// topology::HETEROGENEOUS_FACE -- topology::HETEROGENEOUS_FACE
//***************************************************************************

template <>
struct topology_data<topology::HETEROGENEOUS_FACE>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static const topology::topology_t value = topology::HETEROGENEOUS_FACE;
  static const topology::topology_t base = value;
  static const bool is_valid = true;
  static const topology::rank_t rank = topology::FACE_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true // 3d
                              > spatial_dimension_vector;

};

//***************************************************************************
// topology::HETEROGENEOUS_ELEMENT_2D -- topology::HETEROGENEOUS_ELEMENT_2D
//***************************************************************************

template <>
struct topology_data<topology::HETEROGENEOUS_ELEMENT_2D>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static const topology::topology_t value = topology::HETEROGENEOUS_ELEMENT_2D;
  static const topology::topology_t base = value;
  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true // 2d
                                , false // 3d
                              > spatial_dimension_vector;

};

//***************************************************************************
// topology::HETEROGENEOUS_ELEMENT -- topology::HETEROGENEOUS_ELEMENT
//***************************************************************************

template <>
struct topology_data<topology::HETEROGENEOUS_ELEMENT>
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static const topology::topology_t value = topology::HETEROGENEOUS_ELEMENT;
  static const topology::topology_t base = value;
  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true // 3d
                              > spatial_dimension_vector;

};

//***************************************************************************
// topology::NODE -- topology::NODE_RANK
//***************************************************************************

template <>
struct topology_data<topology::NODE>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::NODE;
  static const topology::topology_t base = topology::NODE;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::NODE_RANK;
  static const topology::rank_t side_rank = topology::INVALID_RANK;
  static const topology::topology_t edge_topology = topology::INVALID_TOPOLOGY;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 0;
  static const unsigned num_nodes = 0;
  static const unsigned num_vertices = 0;
  static const unsigned num_edges = 0;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 0;
  static const unsigned num_positive_permutations = 0;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , true  // 1d
                                , true  // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<> edge_node_ordinals_vector;
  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<> permutation_node_ordinals_vector;

};

//***************************************************************************
// PARTICLE -- topology::ELEMENT_RANK
// one node
//***************************************************************************

template <>
struct topology_data<topology::PARTICLE>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::PARTICLE;
  static const topology::topology_t base = topology::PARTICLE;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::NODE_RANK;
  static const topology::topology_t edge_topology = topology::INVALID_TOPOLOGY;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 1;
  static const unsigned num_nodes = 1;
  static const unsigned num_vertices = 1;
  static const unsigned num_edges = 0;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 1;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , true  // 1d
                                , true  // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<> edge_node_ordinals_vector;
  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0>
                            > permutation_node_ordinals_vector;

};


//***************************************************************************
// topology::LINE -- topology::EDGE_RANK
// 2 or 3 nodes
//
//  o------o------o
//  0      2      1
//
//***************************************************************************

template <>
struct topology_data<topology::LINE_2>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::LINE_2;
  static const topology::topology_t base = topology::LINE_2;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::EDGE_RANK;
  static const topology::rank_t side_rank = topology::NODE_RANK;
  static const topology::topology_t edge_topology = topology::INVALID_TOPOLOGY;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 1;
  static const unsigned num_nodes = 2;
  static const unsigned num_vertices = 2;
  static const unsigned num_edges = 0;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 2;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<> edge_node_ordinals_vector;
  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 0>
                            > permutation_node_ordinals_vector;

};


template <>
struct topology_data<topology::LINE_3>
  : public topology_data<topology::LINE_2>
{
  static const topology::topology_t value = topology::LINE_3;
  static const unsigned num_nodes = 3;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2>
    , boost::mpl::vector_c<unsigned, 1, 0, 2>
                            > permutation_node_ordinals_vector;
};

//***************************************************************************
// topology::LINE 1D -- topology::ELEMENT_RANK
// only defined on 1d problems
//
//  o------o------o
//  0      2      1
//
//***************************************************************************

template <>
struct topology_data<topology::LINE_2_1D>
  : public topology_data<topology::LINE_2>
{
  static const topology::topology_t value = topology::LINE_2_1D;
  static const topology::topology_t base = topology::LINE_2_1D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , true  // 1d
                                , false // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

template <>
struct topology_data<topology::LINE_3_1D>
  : public topology_data<topology::LINE_3>
{
  static const topology::topology_t value = topology::LINE_3_1D;
  static const topology::topology_t base = topology::LINE_2_1D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , true  // 1d
                                , false // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

//***************************************************************************
// topology::BEAM -- topology::ELEMENT_RANK
// 2 or 3 nodes with a single edge
//
//  o------o------o
//  0      2      1
//       Edge 0
//***************************************************************************

template <>
struct topology_data<topology::BEAM_2>
  : public topology_data<topology::LINE_2>
{
  static const topology::topology_t value = topology::BEAM_2;
  static const topology::topology_t base = topology::BEAM_2;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;

  static const bool is_shell = false;
  static const unsigned dimension = 2;
  static const unsigned num_edges = 1;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
                            > edge_node_ordinals_vector;

};

template <>
struct topology_data<topology::BEAM_3>
  : public topology_data<topology::LINE_3>
{
  static const topology::topology_t value = topology::BEAM_3;
  static const topology::topology_t base = topology::BEAM_2;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_3;

  static const bool is_shell = false;
  static const unsigned dimension = 2;
  static const unsigned num_edges = 1;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2>
                            > edge_node_ordinals_vector;

};

//***************************************************************************
// topology::SHELL_LINE -- topology::ELEMENT_RANK
// only defined on 2d problems
// 2 or 3 nodes with a two edge
//
//       Edge 1: (1,0,2)
//
//  o------o------o
//  0      2      1
//
//       Edge 0: (0,1,2)
//***************************************************************************

template <>
struct topology_data<topology::SHELL_LINE_2>
  : public topology_data<topology::LINE_2>
{
  static const topology::topology_t value = topology::SHELL_LINE_2;
  static const topology::topology_t base = topology::SHELL_LINE_2;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;

  static const bool is_shell = true;
  static const unsigned dimension = 2;
  static const unsigned num_edges = 2;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 0>
                            > edge_node_ordinals_vector;

};

template <>
struct topology_data<topology::SHELL_LINE_3>
  : public topology_data<topology::LINE_3>
{
  static const topology::topology_t value = topology::SHELL_LINE_3;
  static const topology::topology_t base = topology::SHELL_LINE_2;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_3;

  static const bool is_shell = true;
  static const unsigned dimension = 2;
  static const unsigned num_edges = 2;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2>
    , boost::mpl::vector_c<unsigned, 1, 0, 2>
                            > edge_node_ordinals_vector;

};

//***************************************************************************
// topology::TRIANGLE -- topology::FACE_RANK
// defined on spatial dimension 3d
// 3, 4, or 6 nodes with 3 edges
/*
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
*/
//***************************************************************************

template <>
struct topology_data<topology::TRI_3>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::TRI_3;
  static const topology::topology_t base = topology::TRI_3;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::FACE_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 2;
  static const unsigned num_nodes = 3;
  static const unsigned num_vertices = 3;
  static const unsigned num_edges = 3;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 6;
  static const unsigned num_positive_permutations = 3;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 0>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2, 0>
    , boost::mpl::vector_c<unsigned, 0, 2, 1>
    , boost::mpl::vector_c<unsigned, 2, 1, 0>
    , boost::mpl::vector_c<unsigned, 1, 0, 2>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::TRI_4>
  : public topology_data<topology::TRI_3>
{
  static const topology::topology_t value = topology::TRI_4;
  static const unsigned num_nodes = 4;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2,  3>
    , boost::mpl::vector_c<unsigned, 2, 0, 1,  3>
    , boost::mpl::vector_c<unsigned, 1, 2, 0,  3>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,  3>
    , boost::mpl::vector_c<unsigned, 2, 1, 0,  3>
    , boost::mpl::vector_c<unsigned, 1, 0, 2,  3>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::TRI_6>
  : public topology_data<topology::TRI_3>
{
  static const topology::topology_t value = topology::TRI_6;
  static const unsigned num_nodes = 6;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  3>
    , boost::mpl::vector_c<unsigned, 1, 2,  4>
    , boost::mpl::vector_c<unsigned, 2, 0,  5>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2,  3, 4, 5>
    , boost::mpl::vector_c<unsigned, 2, 0, 1,  5, 3, 4>
    , boost::mpl::vector_c<unsigned, 1, 2, 0,  4, 5, 3>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,  5, 4, 3>
    , boost::mpl::vector_c<unsigned, 2, 1, 0,  4, 3, 5>
    , boost::mpl::vector_c<unsigned, 1, 0, 2,  3, 5, 4>
                            > permutation_node_ordinals_vector;
};

//***************************************************************************
// topology::TRIANGLE 2D -- topology::ELEMENT_RANK
// defined on spatial dimension 2d
// 3, 4, or 6 nodes with 3 edges
/*
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
*/
//***************************************************************************

template <>
struct topology_data<topology::TRI_3_2D>
  : public topology_data<topology::TRI_3>
{
  static const topology::topology_t value = topology::TRI_3_2D;
  static const topology::topology_t base = topology::TRI_3_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

template <>
struct topology_data<topology::TRI_4_2D>
  : public topology_data<topology::TRI_4>
{
  static const topology::topology_t value = topology::TRI_4_2D;
  static const topology::topology_t base = topology::TRI_3_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

template <>
struct topology_data<topology::TRI_6_2D>
  : public topology_data<topology::TRI_6>
{
  static const topology::topology_t value = topology::TRI_6_2D;
  static const topology::topology_t base = topology::TRI_3_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

//***************************************************************************
// topology::SHELL topology::TRIANGLE -- topology::ELEMENT_RANK
// defined on spatial dimension 3d
// 3, 4, or 6 nodes with 3 edges and 2 faces
/*
//
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2    5 o       o 4   Edge #1
//               /         \
//              /           \
//             /             \
//            o-------o-------o
//           0        3        1
//
//                  Edge #0
//
//   Face #0 (0, 1, 2,   3, 4, 5)
//   Face #1 (0, 2, 1,   5, 4, 3)
//
//
//                    2
//                    o
//                   / \
//                  /   \
//                 /     \
//   Edge #2      /       \     Edge #1
//               /    4    \
//              /     o     \
//             /             \
//            o---------------o
//           0                 1
//
//                  Edge #0
//
//   Face #0 (0, 1, 2,   3)
//   Face #1 (0, 2, 1,   3)
*/
//***************************************************************************

template <>
struct topology_data<topology::SHELL_TRI_3>
  : public topology_data<topology::TRI_3>
{
  static const topology::topology_t value = topology::SHELL_TRI_3;
  static const topology::topology_t base = topology::SHELL_TRI_3;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_3
                                , topology::TRI_3
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2>
    , boost::mpl::vector_c<unsigned, 0, 2, 1>
                            > face_node_ordinals_vector;
};

template <>
struct topology_data<topology::SHELL_TRI_4>
  : public topology_data<topology::TRI_4>
{
  static const topology::topology_t value = topology::SHELL_TRI_4;
  static const topology::topology_t base = topology::SHELL_TRI_3;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_4
                                , topology::TRI_4
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3>
    , boost::mpl::vector_c<unsigned, 0, 2, 1, 3>
                            > face_node_ordinals_vector;
};

template <>
struct topology_data<topology::SHELL_TRI_6>
  : public topology_data<topology::TRI_6>
{
  static const topology::topology_t value = topology::SHELL_TRI_6;
  static const topology::topology_t base = topology::SHELL_TRI_3;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_6
                                , topology::TRI_6
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5>
    , boost::mpl::vector_c<unsigned, 0, 2, 1, 5, 4, 3>
                            > face_node_ordinals_vector;
};

//***************************************************************************
// topology::QUADRILATERAL -- topology::FACE_RANK
// defined on spatial dimension 3d
// 4, 8, or 9 nodes with 4 edges
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//***************************************************************************

template <>
struct topology_data<topology::QUAD_4>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::QUAD_4;
  static const topology::topology_t base = topology::QUAD_4;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::FACE_RANK;
  static const topology::rank_t side_rank = topology::EDGE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 2;
  static const unsigned num_nodes = 4;
  static const unsigned num_vertices = 4;
  static const unsigned num_edges = 4;
  static const unsigned num_faces = 0;
  static const unsigned num_permutations = 8;
  static const unsigned num_positive_permutations = 4;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<topology::topology_t> face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 3>
    , boost::mpl::vector_c<unsigned, 3, 0>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<> face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3>
    , boost::mpl::vector_c<unsigned, 3, 0, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 3, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2, 3, 0>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1>
    , boost::mpl::vector_c<unsigned, 3, 2, 1, 0>
    , boost::mpl::vector_c<unsigned, 2, 1, 0, 3>
    , boost::mpl::vector_c<unsigned, 1, 0, 3, 2>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::QUAD_8>
  : public topology_data<topology::QUAD_4>
{
  static const topology::topology_t value = topology::QUAD_8;
  static const unsigned num_nodes = 8;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  4>
    , boost::mpl::vector_c<unsigned, 1, 2,  5>
    , boost::mpl::vector_c<unsigned, 2, 3,  6>
    , boost::mpl::vector_c<unsigned, 3, 0,  7>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7>
    , boost::mpl::vector_c<unsigned, 3, 0, 1, 2,  7, 4, 5, 6>
    , boost::mpl::vector_c<unsigned, 2, 3, 0, 1,  6, 7, 4, 5>
    , boost::mpl::vector_c<unsigned, 1, 2, 3, 0,  5, 6, 7, 4>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  7, 6, 5, 4>
    , boost::mpl::vector_c<unsigned, 3, 2, 1, 0,  6, 5, 4, 7>
    , boost::mpl::vector_c<unsigned, 2, 1, 0, 3,  5, 4, 7, 6>
    , boost::mpl::vector_c<unsigned, 1, 0, 3, 2,  4, 7, 6, 5>
                            > permutation_node_ordinals_vector;
};

template <>
struct topology_data<topology::QUAD_9>
  : public topology_data<topology::QUAD_8>
{
  static const topology::topology_t value = topology::QUAD_9;
  static const unsigned num_nodes = 9;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7,  8>
    , boost::mpl::vector_c<unsigned, 3, 0, 1, 2,  7, 4, 5, 6,  8>
    , boost::mpl::vector_c<unsigned, 2, 3, 0, 1,  6, 7, 4, 5,  8>
    , boost::mpl::vector_c<unsigned, 1, 2, 3, 0,  5, 6, 7, 4,  8>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  7, 6, 5, 4,  8>
    , boost::mpl::vector_c<unsigned, 3, 2, 1, 0,  6, 5, 4, 7,  8>
    , boost::mpl::vector_c<unsigned, 2, 1, 0, 3,  5, 4, 7, 6,  8>
    , boost::mpl::vector_c<unsigned, 1, 0, 3, 2,  4, 7, 6, 5,  8>
                            > permutation_node_ordinals_vector;
};

//***************************************************************************
// topology::QUADRILATERAL 2D -- topology::ELEMENT_RANK
// defined on spatial dimension 2d
// 4, 8, or 9 nodes with 4 edges
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//***************************************************************************

template <>
struct topology_data<topology::QUAD_4_2D>
  : public topology_data<topology::QUAD_4>
{
  static const topology::topology_t value = topology::QUAD_4_2D;
  static const topology::topology_t base = topology::QUAD_4_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

template <>
struct topology_data<topology::QUAD_8_2D>
  : public topology_data<topology::QUAD_8>
{
  static const topology::topology_t value = topology::QUAD_8_2D;
  static const topology::topology_t base = topology::QUAD_4_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

template <>
struct topology_data<topology::QUAD_9_2D>
  : public topology_data<topology::QUAD_9>
{
  static const topology::topology_t value = topology::QUAD_9_2D;
  static const topology::topology_t base = topology::QUAD_4_2D;

  static const topology::rank_t rank = topology::ELEMENT_RANK;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , true  // 2d
                                , false // 3d
                              > spatial_dimension_vector;
};

//***************************************************************************
//  topology::SHELL topology::QUADRILATERAL -- topology::ELEMENT_RANK
// defined on spatial dimension 3d
// 4, 8, or 9 nodes with 4 edges and 2 faces
//
//                 Edge #2
//
//            3        6        2
//             o-------o-------o
//             |               |
//             |               |
//             |       8       |
//  Edge #3  7 o       o       o 5  Edge #1
//             |               |
//             |               |
//             |               |
//             o-------o-------o
//            0        4        1
//
//                  Edge #0
//
//  Face #0 (0, 1, 2, 3,   4, 5, 6, 7,   8)
//  Face #1 (0, 3, 2, 1,   7, 6, 5, 4,   8)
//
//***************************************************************************

template <>
struct topology_data<topology::SHELL_QUAD_4>
  : public topology_data<topology::QUAD_4>
{
  static const topology::topology_t value = topology::SHELL_QUAD_4;
  static const topology::topology_t base = topology::SHELL_QUAD_4;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_4
                                , topology::QUAD_4
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1>
                            > face_node_ordinals_vector;
};

template <>
struct topology_data<topology::SHELL_QUAD_8>
  : public topology_data<topology::QUAD_8>
{
  static const topology::topology_t value = topology::SHELL_QUAD_8;
  static const topology::topology_t base = topology::SHELL_QUAD_4;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_8
                                , topology::QUAD_8
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  7, 6, 5, 4>
                            > face_node_ordinals_vector;
};

template <>
struct topology_data<topology::SHELL_QUAD_9>
  : public topology_data<topology::QUAD_9>
{
  static const topology::topology_t value = topology::SHELL_QUAD_9;
  static const topology::topology_t base = topology::SHELL_QUAD_4;

  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const bool is_shell = true;

  static const unsigned dimension = 3;
  static const unsigned num_faces = 2;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_9
                                , topology::QUAD_9
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7,  8>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  7, 6, 5, 4,  8>
                            > face_node_ordinals_vector;
};

//***************************************************************************
// topology::TETRAHEDRON
//***************************************************************************

template <>
struct topology_data<topology::TET_4>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::TET_4;
  static const topology::topology_t base = topology::TET_4;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = true;
  static const bool is_shell = false;
  static const unsigned dimension = 3;
  static const unsigned num_nodes = 4;
  static const unsigned num_vertices = 4;
  static const unsigned num_edges = 6;
  static const unsigned num_faces = 4;
  static const unsigned num_permutations = 1;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_3
                                , topology::TRI_3
                                , topology::TRI_3
                                , topology::TRI_3
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 0>
    , boost::mpl::vector_c<unsigned, 0, 3>
    , boost::mpl::vector_c<unsigned, 1, 3>
    , boost::mpl::vector_c<unsigned, 2, 3>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 3>
    , boost::mpl::vector_c<unsigned, 1, 2, 3>
    , boost::mpl::vector_c<unsigned, 0, 3, 2>
    , boost::mpl::vector_c<unsigned, 0, 2, 1>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::TET_8>
  : public topology_data<topology::TET_4>
{
  static const topology::topology_t value = topology::TET_8;
  static const unsigned num_nodes = 8;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_4
                                , topology::TRI_4
                                , topology::TRI_4
                                , topology::TRI_4
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 3,  4>
    , boost::mpl::vector_c<unsigned, 1, 2, 3,  5>
    , boost::mpl::vector_c<unsigned, 0, 3, 2,  7>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,  6>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::TET_10>
  : public topology_data<topology::TET_4>
{
  static const topology::topology_t value = topology::TET_10;
  static const unsigned num_nodes = 10;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  4>
    , boost::mpl::vector_c<unsigned, 1, 2,  5>
    , boost::mpl::vector_c<unsigned, 2, 0,  6>
    , boost::mpl::vector_c<unsigned, 0, 3,  7>
    , boost::mpl::vector_c<unsigned, 1, 3,  8>
    , boost::mpl::vector_c<unsigned, 2, 3,  9>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 3,  4, 8, 7>
    , boost::mpl::vector_c<unsigned, 1, 2, 3,  5, 9, 8>
    , boost::mpl::vector_c<unsigned, 0, 3, 2,  7, 9, 6>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,  6, 5, 4>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7, 8, 9>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::TET_11>
  : public topology_data<topology::TET_10>
{
  static const topology::topology_t value = topology::TET_11;
  static const unsigned num_nodes = 11;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3,  4, 5, 6, 7, 8, 9,  10>
                            > permutation_node_ordinals_vector;
};

//***************************************************************************
// topology::PYRAMID
//***************************************************************************
template <>
struct topology_data<topology::PYRAMID_5>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::PYRAMID_5;
  static const topology::topology_t base = topology::PYRAMID_5;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 3;
  static const unsigned num_nodes = 5;
  static const unsigned num_vertices = 5;
  static const unsigned num_edges = 8;
  static const unsigned num_faces = 5;
  static const unsigned num_permutations = 1;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_3
                                , topology::TRI_3
                                , topology::TRI_3
                                , topology::TRI_3
                                , topology::QUAD_4
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 3>
    , boost::mpl::vector_c<unsigned, 3, 0>
    , boost::mpl::vector_c<unsigned, 0, 4>
    , boost::mpl::vector_c<unsigned, 1, 4>
    , boost::mpl::vector_c<unsigned, 2, 4>
    , boost::mpl::vector_c<unsigned, 3, 4>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4>
    , boost::mpl::vector_c<unsigned, 1, 2, 4>
    , boost::mpl::vector_c<unsigned, 2, 3, 4>
    , boost::mpl::vector_c<unsigned, 3, 0, 4>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::PYRAMID_13>
  : public topology_data<topology::PYRAMID_5>
{
  static const topology::topology_t value = topology::PYRAMID_13;
  static const unsigned num_nodes = 13;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::QUAD_8
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  5>
    , boost::mpl::vector_c<unsigned, 1, 2,  6>
    , boost::mpl::vector_c<unsigned, 2, 3,  7>
    , boost::mpl::vector_c<unsigned, 3, 0,  8>
    , boost::mpl::vector_c<unsigned, 0, 4,  9>
    , boost::mpl::vector_c<unsigned, 1, 4,  10>
    , boost::mpl::vector_c<unsigned, 2, 4,  11>
    , boost::mpl::vector_c<unsigned, 3, 4,  12>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4,  5, 10, 9>
    , boost::mpl::vector_c<unsigned, 1, 2, 4,  6, 11, 10>
    , boost::mpl::vector_c<unsigned, 2, 3, 4,  7, 12, 11>
    , boost::mpl::vector_c<unsigned, 3, 0, 4,  8, 9,  12>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  8, 7, 6, 5>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11, 12>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::PYRAMID_14>
  : public topology_data<topology::PYRAMID_13>
{
  static const topology::topology_t value = topology::PYRAMID_14;
  static const unsigned num_nodes = 14;


  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::TRI_6
                                , topology::QUAD_9
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4,  5, 10, 9>
    , boost::mpl::vector_c<unsigned, 1, 2, 4,  6, 11, 10>
    , boost::mpl::vector_c<unsigned, 2, 3, 4,  7, 12, 11>
    , boost::mpl::vector_c<unsigned, 3, 0, 4,  8, 9,  12>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  8, 7, 6, 5,  13>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4,  5, 6, 7, 8, 9, 10, 11, 12,  13>
                            > permutation_node_ordinals_vector;

};

//***************************************************************************
// topology::WEDGE
//***************************************************************************
template <>
struct topology_data<topology::WEDGE_6>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::WEDGE_6;
  static const topology::topology_t base = topology::WEDGE_6;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = false;
  static const bool is_shell = false;
  static const unsigned dimension = 3;
  static const unsigned num_nodes = 6;
  static const unsigned num_vertices = 6;
  static const unsigned num_edges = 9;
  static const unsigned num_faces = 5;
  static const unsigned num_permutations = 1;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::TRI_3
                                , topology::TRI_3
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 0>
    , boost::mpl::vector_c<unsigned, 3, 4>
    , boost::mpl::vector_c<unsigned, 4, 5>
    , boost::mpl::vector_c<unsigned, 5, 3>
    , boost::mpl::vector_c<unsigned, 0, 3>
    , boost::mpl::vector_c<unsigned, 1, 4>
    , boost::mpl::vector_c<unsigned, 2, 5>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4, 3>
    , boost::mpl::vector_c<unsigned, 1, 2, 5, 4>
    , boost::mpl::vector_c<unsigned, 0, 3, 5, 2>
    , boost::mpl::vector_c<unsigned, 0, 2, 1>
    , boost::mpl::vector_c<unsigned, 3, 4, 5>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::WEDGE_15>
  : public topology_data<topology::WEDGE_6>
{
  static const topology::topology_t value = topology::WEDGE_15;
  static const unsigned num_nodes = 15;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::TRI_6
                                , topology::TRI_6
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  6>
    , boost::mpl::vector_c<unsigned, 1, 2,  7>
    , boost::mpl::vector_c<unsigned, 2, 0,  8>
    , boost::mpl::vector_c<unsigned, 3, 4,  12>
    , boost::mpl::vector_c<unsigned, 4, 5,  13>
    , boost::mpl::vector_c<unsigned, 5, 3,  14>
    , boost::mpl::vector_c<unsigned, 0, 3,  9>
    , boost::mpl::vector_c<unsigned, 1, 4,  10>
    , boost::mpl::vector_c<unsigned, 2, 5,  11>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4, 3,  6, 10, 12,  9>
    , boost::mpl::vector_c<unsigned, 1, 2, 5, 4,  7, 11, 13, 10>
    , boost::mpl::vector_c<unsigned, 0, 3, 5, 2,  9, 14, 11,  8>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,   8,  7, 6>
    , boost::mpl::vector_c<unsigned, 3, 4, 5,  12, 13, 14>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5, 6,  7, 8, 9, 10, 11, 12, 13, 14>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::WEDGE_18>
  : public topology_data<topology::WEDGE_15>
{
  static const topology::topology_t value = topology::WEDGE_18;
  static const unsigned num_nodes = 18;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::TRI_6
                                , topology::TRI_6
                              > face_topology_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 4, 3,  6, 10, 12,  9,  15>
    , boost::mpl::vector_c<unsigned, 1, 2, 5, 4,  7, 11, 13, 10,  16>
    , boost::mpl::vector_c<unsigned, 0, 3, 5, 2,  9, 14, 11,  8,  17>
    , boost::mpl::vector_c<unsigned, 0, 2, 1,   8,  7, 6>
    , boost::mpl::vector_c<unsigned, 3, 4, 5,  12, 13, 14>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5, 6,  7, 8, 9, 10, 11, 12, 13, 14,  15, 16, 17>
                            > permutation_node_ordinals_vector;

};

//***************************************************************************
// topology::HEXAHEDRON -- topology::ELEMENT_RANK
// defined on 3d
/*
//   Linear 8-Node Hexahedron node locations.
//
//          7                    6
//           o------------------o
//          /|                 /|
//         / |                / |
//        /  |               /  |
//       /   |              /   |
//      /    |             /    |
//     /     |            /     |
//  4 /      |         5 /      |
//   o------------------o       |
//   |       |          |       |
//   |     3 o----------|-------o 2
//   |      /           |      /
//   |     /            |     /
//   |    /             |    /
//   |   /              |   /
//   |  /               |  /
//   | /                | /
//   |/                 |/
//   o------------------o
//  0                    1
//
//
//   Quadratic 20-Node Hexahedron node locations:
//
//           7         18         6
//            o--------o---------o
//           /|                 /|
//          / |                / |
//         /  |               /  |
//      19o   |            17o   |
//       /  15o             /    o14
//      /     |            /     |
//   4 /      | 16        /      |
//    o---------o--------o 5     |
//    |       |       10 |       |
//    |     3 o-------o--|-------o 2
//    |      /           |      /
//    |     /            |     /
//  12o    /             o13  /
//    |   o11            |   o9
//    |  /               |  /
//    | /                | /
//    |/                 |/
//    o---------o--------o
//   0          8         1
//
//
//   Quadratic 27-Node Hexahedron additional node locations:
//
//
//            x--------x---------x
//           /|                 /|
//          / |                / |
//         /  |   22          /  |
//        x   |    o         x   |
//       /    x       o26   /    x     Node #20 is at centroid of element
//      /     |            /     |
//     /      |           /      |     "QUAD_9" beginning with nodes
//    x---------x--------x       |      0,1,5,4 has node 25 at center....
//    | 23o   |          |   o24 |
//    |       x-------x--|-------x
//    |      /           |      /
//    |     /  25        |     /
//    x    /    o        x    /
//    |   x        o21   |   x
//    |  /               |  /
//    | /                | /
//    |/                 |/
//    x---------x--------x
*/
//***************************************************************************
template <>
struct topology_data<topology::HEX_8>
{
  typedef topology::topology_t value_type;
  static const topology::topology_t value = topology::HEX_8;
  static const topology::topology_t base = topology::HEX_8;

  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::FACE_RANK;
  static const topology::topology_t edge_topology = topology::LINE_2;
  static const bool has_homogeneous_faces = true;
  static const bool is_shell = false;
  static const unsigned dimension = 3;
  static const unsigned num_nodes = 8;
  static const unsigned num_vertices = 8;
  static const unsigned num_edges = 12;
  static const unsigned num_faces = 6;
  static const unsigned num_permutations = 1;
  static const unsigned num_positive_permutations = 1;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , false // 1d
                                , false // 2d
                                , true  // 3d
                              > spatial_dimension_vector;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::QUAD_4
                                , topology::QUAD_4
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1>
    , boost::mpl::vector_c<unsigned, 1, 2>
    , boost::mpl::vector_c<unsigned, 2, 3>
    , boost::mpl::vector_c<unsigned, 3, 0>
    , boost::mpl::vector_c<unsigned, 4, 5>
    , boost::mpl::vector_c<unsigned, 5, 6>
    , boost::mpl::vector_c<unsigned, 6, 7>
    , boost::mpl::vector_c<unsigned, 7, 4>
    , boost::mpl::vector_c<unsigned, 0, 4>
    , boost::mpl::vector_c<unsigned, 1, 5>
    , boost::mpl::vector_c<unsigned, 2, 6>
    , boost::mpl::vector_c<unsigned, 3, 7>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 5, 4>
    , boost::mpl::vector_c<unsigned, 1, 2, 6, 5>
    , boost::mpl::vector_c<unsigned, 2, 3, 7, 6>
    , boost::mpl::vector_c<unsigned, 0, 4, 7, 3>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1>
    , boost::mpl::vector_c<unsigned, 4, 5, 6, 7>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5, 6, 7>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::HEX_20>
  : public topology_data<topology::HEX_8>
{
  static const topology::topology_t value = topology::HEX_20;
  static const unsigned num_nodes = 20;

  static const topology::topology_t edge_topology = topology::LINE_3;

  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::QUAD_8
                                , topology::QUAD_8
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1,  8>
    , boost::mpl::vector_c<unsigned, 1, 2,  9>
    , boost::mpl::vector_c<unsigned, 2, 3,  10>
    , boost::mpl::vector_c<unsigned, 3, 0,  11>
    , boost::mpl::vector_c<unsigned, 4, 5,  16>
    , boost::mpl::vector_c<unsigned, 5, 6,  17>
    , boost::mpl::vector_c<unsigned, 6, 7,  18>
    , boost::mpl::vector_c<unsigned, 7, 4,  19>
    , boost::mpl::vector_c<unsigned, 0, 4,  12>
    , boost::mpl::vector_c<unsigned, 1, 5,  13>
    , boost::mpl::vector_c<unsigned, 2, 6,  14>
    , boost::mpl::vector_c<unsigned, 3, 7,  15>
                            > edge_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 5, 4,   8, 13, 16, 12>
    , boost::mpl::vector_c<unsigned, 1, 2, 6, 5,   9, 14, 17, 13>
    , boost::mpl::vector_c<unsigned, 2, 3, 7, 6,  10, 15, 18, 14>
    , boost::mpl::vector_c<unsigned, 0, 4, 7, 3,  12, 19, 15, 11>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  11, 10,  9,  8>
    , boost::mpl::vector_c<unsigned, 4, 5, 6, 7,  16, 17, 18, 19>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19>
                            > permutation_node_ordinals_vector;

};

template <>
struct topology_data<topology::HEX_27>
  : public topology_data<topology::HEX_20>
{
  static const topology::topology_t value = topology::HEX_27;
  static const unsigned num_nodes = 27;


  typedef boost::mpl::vector_c<   topology::topology_t
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::QUAD_9
                                , topology::QUAD_9
                              > face_topology_vector;


  typedef boost::mpl::vector<
      boost::mpl::vector_c<unsigned, 0, 1, 5, 4,   8, 13, 16, 12,  25>
    , boost::mpl::vector_c<unsigned, 1, 2, 6, 5,   9, 14, 17, 13,  24>
    , boost::mpl::vector_c<unsigned, 2, 3, 7, 6,  10, 15, 18, 14,  26>
    , boost::mpl::vector_c<unsigned, 0, 4, 7, 3,  12, 19, 15, 11,  23>
    , boost::mpl::vector_c<unsigned, 0, 3, 2, 1,  11, 10,  9,  8,  21>
    , boost::mpl::vector_c<unsigned, 4, 5, 6, 7,  16, 17, 18, 19,  22>
                            > face_node_ordinals_vector;

  typedef boost::mpl::vector<
      boost::mpl::vector27_c<unsigned, 0, 1, 2, 3, 4, 5, 6, 7,  8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  20, 21, 22, 23, 24, 25, 26>
                            > permutation_node_ordinals_vector;

};

//***************************************************************************
// topology::SUPERELEMENT -- topology::SUPERELEMENT
//***************************************************************************

template <topology::topology_t Topology>
struct topology_data<Topology, typename boost::enable_if_c< (Topology > topology::SUPERELEMENT_START) >::type >
  : public topology_data<topology::INVALID_TOPOLOGY>
{
  static const topology::topology_t value = Topology;
  static const topology::topology_t base = value;
  static const bool is_valid = true;
  static const topology::rank_t rank = topology::ELEMENT_RANK;
  static const topology::rank_t side_rank = topology::INVALID_RANK;
  static const unsigned num_nodes = Topology - topology::SUPERELEMENT_START;

  typedef boost::mpl::vector_c<   bool
                                , false // 0d
                                , true // 1d
                                , true // 2d
                                , true // 3d
                              > spatial_dimension_vector;

};

}} // namespace stk::topology_detail

#endif //STKTOPOLOGY_DETAIL_TOPOLOGY_DATA_HPP
