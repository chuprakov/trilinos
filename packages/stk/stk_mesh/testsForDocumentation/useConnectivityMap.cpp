#include <gtest/gtest.h>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ConnectivityMap.hpp>

namespace {

TEST(stkMeshHowTo, useConnectivityMap)
{
    const unsigned spatialDimension = 2;
    stk::mesh::MetaData metaData(spatialDimension, stk::mesh::entity_rank_names());
    metaData.commit();

    const stk::mesh::ConnectivityMap& minimal_connectivity = stk::mesh::ConnectivityMap::minimal_upward_connectivity_map();

    //The minimal_upward_connectivity_map should disable node->edge and edge->elem connections
    //but still allow node->elem connections:
    EXPECT_TRUE(minimal_connectivity.valid(stk::topology::ELEM_RANK, stk::topology::NODE_RANK));
    EXPECT_TRUE(minimal_connectivity.valid(stk::topology::NODE_RANK, stk::topology::ELEM_RANK));
    EXPECT_FALSE(minimal_connectivity.valid(stk::topology::NODE_RANK, stk::topology::EDGE_RANK));
    EXPECT_FALSE(minimal_connectivity.valid(stk::topology::EDGE_RANK, stk::topology::ELEM_RANK));

    bool add_fmwk_data = false;
    stk::mesh::BulkData mesh(metaData, MPI_COMM_WORLD, add_fmwk_data, &minimal_connectivity);
    mesh.modification_begin();

    //set up 1 element (3-node triangle) that is fully downward connected to edges and nodes
    stk::mesh::EntityId elemId = 1;
    stk::mesh::EntityId elemNodeIds[] = {1, 2, 3};
    stk::mesh::EntityId elemEdgeIds[] = {6, 7, 8};
    stk::mesh::Entity elemNodes[3];
    stk::mesh::Entity elemEdges[3];
    stk::mesh::Entity elem = mesh.declare_entity(stk::topology::ELEM_RANK, elemId);
    elemNodes[0] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[0]);
    elemNodes[1] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[1]);
    elemNodes[2] = mesh.declare_entity(stk::topology::NODE_RANK, elemNodeIds[2]);

    elemEdges[0] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[0]);
    elemEdges[1] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[1]);
    elemEdges[2] = mesh.declare_entity(stk::topology::EDGE_RANK, elemEdgeIds[2]);

    //downward element -> node connectivity
    mesh.declare_relation(elem, elemNodes[0], 0);
    mesh.declare_relation(elem, elemNodes[1], 1);
    mesh.declare_relation(elem, elemNodes[2], 2);

    //downward element -> edge connectivity
    mesh.declare_relation(elem, elemEdges[0], 0);
    mesh.declare_relation(elem, elemEdges[1], 1);
    mesh.declare_relation(elem, elemEdges[2], 2);

    //downward edge -> node connectivity
    mesh.declare_relation(elemEdges[0], elemNodes[0], 0); mesh.declare_relation(elemEdges[0], elemNodes[1], 1);
    mesh.declare_relation(elemEdges[1], elemNodes[1], 0); mesh.declare_relation(elemEdges[1], elemNodes[2], 1);
    mesh.declare_relation(elemEdges[2], elemNodes[2], 0); mesh.declare_relation(elemEdges[2], elemNodes[0], 1);
    mesh.modification_end();

    //now test upward connectivity which stk-mesh generally creates by default
    unsigned expectedNumElemsPerEdge = 0;//this would be 1 by default
    unsigned expectedNumEdgesPerNode = 0;//this would be 2 by default
    EXPECT_EQ(expectedNumElemsPerEdge, mesh.num_elements(elemEdges[0]));
    EXPECT_EQ(expectedNumEdgesPerNode, mesh.num_edges(elemNodes[0]));

    //even though stk-mesh didn't store node->edge connectivity (because we disabled it),
    //we can still get it if we need it:
    stk::mesh::EntityVector connected_edges;
    unsigned numEdges = stk::mesh::get_connectivity(mesh, elemNodes[0], stk::topology::EDGE_RANK, connected_edges);
    EXPECT_EQ(2u, numEdges);
    EXPECT_EQ(numEdges, connected_edges.size());
}

}

