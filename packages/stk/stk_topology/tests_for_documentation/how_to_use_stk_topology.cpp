#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <vector>

namespace {

void verifyPermutationsForTriangle(stk::topology triangular_shell, unsigned* triangle_1_node_ids, unsigned* gold_triangle_1_permutations)
{
    ASSERT_TRUE(stk::topology::SHELL_TRIANGLE_3 == triangular_shell);
    unsigned triangle_1_permutation[3];
    for (unsigned i=0;i<triangular_shell.num_permutations();i++)
    {
        triangular_shell.permutation_nodes(triangle_1_node_ids, i, triangle_1_permutation);
        EXPECT_TRUE(gold_triangle_1_permutations[3*i+0] == triangle_1_permutation[0] &&
                    gold_triangle_1_permutations[3*i+1] == triangle_1_permutation[1] &&
                    gold_triangle_1_permutations[3*i+2] == triangle_1_permutation[2]);
    }
}

//Lexicographical
TEST(stk_topology_understanding, lexicographical_smallest_permutation)
{
    {
        unsigned triangle_node_ids[3] = {10, 8, 12};

        stk::topology triangular_shell = stk::topology::SHELL_TRIANGLE_3;

        unsigned gold_triangle_permutations[18]= {
                10, 8, 12,
                12, 10, 8,
                8, 12, 10, // lexicographical smallest permutation by node ids if considering only positive permutations
                10, 12, 8,
                12, 8, 10,
                8, 10, 12  // lexicographical smallest permutation by node ids if considering all permutations
        };

        verifyPermutationsForTriangle(triangular_shell, triangle_node_ids, gold_triangle_permutations);

        bool usePositivePermutationsOnly = false;
        unsigned permutation_index = triangular_shell.lexicographical_smallest_permutation(triangle_node_ids, usePositivePermutationsOnly);
        unsigned gold_lexicographical_smallest_permutation_index = 5;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);

        usePositivePermutationsOnly = true;
        permutation_index = triangular_shell.lexicographical_smallest_permutation(triangle_node_ids, usePositivePermutationsOnly);
        gold_lexicographical_smallest_permutation_index = 2;
        // driven by vertices, NOT mid-edge nodes
        EXPECT_EQ(gold_lexicographical_smallest_permutation_index, permutation_index);
    }
}

//SubTopology
TEST(stk_topology_understanding, sub_topology)
{
    stk::topology hex20 = stk::topology::HEX_20;
    unsigned hex20Nodes[20] = {
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9, 10, 11,
            12, 13, 14, 15,
            16, 17, 18, 19
    };

    unsigned numFaces = hex20.num_sub_topology(stk::topology::FACE_RANK);
    EXPECT_EQ(6u, numFaces);

    unsigned faceIndex=2;
    stk::topology top = hex20.sub_topology(stk::topology::FACE_RANK, faceIndex);
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, top);

    unsigned nodeIdsFace[8];
    hex20.sub_topology_nodes(hex20Nodes, stk::topology::FACE_RANK, faceIndex, nodeIdsFace);

    unsigned goldIdsFace[8] = { 2, 3, 7, 6, 10, 15, 18, 14 };
    for (unsigned i=0;i<hex20.face_topology(faceIndex).num_nodes();i++)
    {
        EXPECT_EQ(goldIdsFace[i], nodeIdsFace[i]);
    }
}

//Sides
TEST(stk_topology_understanding, sides)
{
    stk::topology hex20 = stk::topology::HEX_20;
    EXPECT_EQ(6u, hex20.num_sides());

    stk::topology quad8 = stk::topology::SHELL_QUADRILATERAL_8;
    EXPECT_EQ(2u, quad8.num_sides());

    stk::topology wedge = stk::topology::WEDGE_15;
    EXPECT_EQ(5u, wedge.num_sides());
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(0));
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(1));
    EXPECT_EQ(stk::topology::QUADRILATERAL_8, wedge.side_topology(2));
    EXPECT_EQ(stk::topology::TRIANGLE_6, wedge.side_topology(3));
    EXPECT_EQ(stk::topology::TRIANGLE_6, wedge.side_topology(4));

}

//Superelements
TEST(stk_topology_understanding, superelements)
{
    unsigned eightNodes=8;
    stk::topology validSuperElement = stk::create_superelement_topology(eightNodes);
    EXPECT_TRUE(validSuperElement.is_superelement());
    EXPECT_TRUE(stk::topology::ELEMENT_RANK == validSuperElement.rank());
    EXPECT_EQ(eightNodes, validSuperElement.num_nodes());
    EXPECT_EQ(0u, validSuperElement.num_edges());
    EXPECT_EQ(0u, validSuperElement.num_faces());
    EXPECT_EQ(0u, validSuperElement.num_permutations());
    EXPECT_EQ(0u, validSuperElement.num_sides());
    EXPECT_EQ(0u, validSuperElement.dimension());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.face_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.edge_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, validSuperElement.base());
    EXPECT_FALSE(validSuperElement.has_homogeneous_faces());
    EXPECT_FALSE(validSuperElement.is_shell());

    unsigned zeroNodes=0;
    stk::topology invalidSuperElement = stk::create_superelement_topology(zeroNodes);
    EXPECT_FALSE(invalidSuperElement.is_superelement());
    EXPECT_TRUE(stk::topology::INVALID_RANK == invalidSuperElement.rank());
    EXPECT_EQ(zeroNodes, invalidSuperElement.num_nodes());
    EXPECT_EQ(0u, invalidSuperElement.num_edges());
    EXPECT_EQ(0u, invalidSuperElement.num_faces());
    EXPECT_EQ(0u, invalidSuperElement.num_permutations());
    EXPECT_EQ(0u, invalidSuperElement.num_sides());
    EXPECT_EQ(0u, invalidSuperElement.dimension());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.face_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.edge_topology());
    EXPECT_EQ(stk::topology::INVALID_TOPOLOGY, invalidSuperElement.base());
    EXPECT_FALSE(invalidSuperElement.has_homogeneous_faces());
    EXPECT_FALSE(invalidSuperElement.is_shell());
}
//Done

}

