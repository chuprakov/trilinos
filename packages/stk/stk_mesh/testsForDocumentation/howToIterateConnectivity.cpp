#include <gtest/gtest.h>                // for AssertHelper, EXPECT_NE, etc
#include <stddef.h>                     // for size_t
#include <limits>                       // for numeric_limits
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_topology/topology.hpp>    // for topology, etc
#include <string>                       // for string
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/CoordinateSystems.hpp"  // for Cartesian
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, FieldBase
#include "stk_mesh/base/Types.hpp"      // for BucketVector
namespace stk { namespace mesh { class Part; } }

namespace
{
//-BEGIN
TEST(StkMeshHowTo, iterateConnectivityThroughBulkData)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // Generate a mesh of hexes with a sideset
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    const stk::mesh::BucketVector &elementBuckets =
        stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);

    typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinatesField_t;
    CoordinatesField_t const & coord_field =
          *dynamic_cast<CoordinatesField_t const *>(stkMeshMetaData.coordinate_field());

    const unsigned nodesPerHex = 8;
    const unsigned spatialDim = 3;
    unsigned count = 0;
    double elementNodeCoords[nodesPerHex][spatialDim];
    for (size_t bucketIndex = 0; bucketIndex < elementBuckets.size(); ++bucketIndex)
    {
        stk::mesh::Bucket &elemBucket = *elementBuckets[bucketIndex];
        for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex)
        {
            stk::mesh::Entity elem = elemBucket[elemIndex];
            unsigned numNodes = stkMeshBulkData.num_nodes(elem);
            EXPECT_EQ(numNodes, nodesPerHex);
            stk::mesh::Entity const* nodes = stkMeshBulkData.begin_nodes(elem);
            for (unsigned inode = 0; inode < numNodes; ++inode)
            {
              double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
              elementNodeCoords[inode][0] = coords[0];
              elementNodeCoords[inode][1] = coords[1];
              elementNodeCoords[inode][2] = coords[2];
              EXPECT_NE(elementNodeCoords[inode][0], std::numeric_limits<double>::max());
              EXPECT_NE(elementNodeCoords[inode][1], std::numeric_limits<double>::max());
              EXPECT_NE(elementNodeCoords[inode][2], std::numeric_limits<double>::max());
              ++count;
            }
        }
    }
    EXPECT_GE(count, 1u);
}

TEST(StkMeshHowTo, iterateConnectivityThroughBuckets)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    stk::io::StkMeshIoBroker stkMeshIoBroker(communicator);
    // Generate a mesh of hexes with a sideset
    const std::string generatedMeshSpecification = "generated:2x2x2|sideset:X";
    stkMeshIoBroker.add_mesh_database(generatedMeshSpecification, stk::io::READ_MESH);
    stkMeshIoBroker.create_input_mesh();
    stkMeshIoBroker.populate_bulk_data();

    stk::mesh::MetaData &stkMeshMetaData = stkMeshIoBroker.meta_data();
    stk::mesh::BulkData &stkMeshBulkData = stkMeshIoBroker.bulk_data();
    const stk::mesh::BucketVector &elementBuckets =
        stkMeshBulkData.buckets(stk::topology::ELEMENT_RANK);

    typedef stk::mesh::Field<double, stk::mesh::Cartesian> CoordinatesField_t;
    CoordinatesField_t const & coord_field =
          *dynamic_cast<CoordinatesField_t const *>(stkMeshMetaData.coordinate_field());

    const unsigned nodesPerHex = 8;
    const unsigned spatialDim = 3;
    unsigned count = 0;
    double elementNodeCoords[nodesPerHex][spatialDim];
    for (size_t bucketIndex = 0; bucketIndex < elementBuckets.size(); ++bucketIndex)
    {
        stk::mesh::Bucket &elemBucket = *elementBuckets[bucketIndex];
        for (size_t elemIndex = 0; elemIndex < elemBucket.size(); ++elemIndex)
        {
            unsigned numNodes = elemBucket.num_nodes(elemIndex);
            EXPECT_EQ(numNodes, nodesPerHex);
            stk::mesh::Entity const* nodes = elemBucket.begin_nodes(elemIndex);
            for (unsigned inode = 0; inode < numNodes; ++inode)
            {
              double *coords = stk::mesh::field_data(coord_field, nodes[inode]);
              elementNodeCoords[inode][0] = coords[0];
              elementNodeCoords[inode][1] = coords[1];
              elementNodeCoords[inode][2] = coords[2];
              EXPECT_NE(elementNodeCoords[inode][0], std::numeric_limits<double>::max());
              EXPECT_NE(elementNodeCoords[inode][1], std::numeric_limits<double>::max());
              EXPECT_NE(elementNodeCoords[inode][2], std::numeric_limits<double>::max());
              ++count;
            }
        }
    }
    EXPECT_GE(count, 1u);
}
//-END
}
