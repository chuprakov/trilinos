/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <vector>                       // for vector, etc

#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>

#include <stdio.h> // getline

//====================
extern int gl_argc;
extern char** gl_argv;

inline std::string getOption(const std::string& option, const std::string defaultString="no")
{
    std::string returnValue = defaultString;
    if ( gl_argv != 0 )
    {
        for (int i=0;i<gl_argc;i++)
        {
            std::string input_argv(gl_argv[i]);
            if ( option == input_argv )
            {
                if ( (i+1) < gl_argc )
                {
                    returnValue = std::string(gl_argv[i+1]);
                }
                break;
            }
        }
    }
    return returnValue;
}

namespace
{

class BulkDataTester : public stk::mesh::BulkData
{
public:
    BulkDataTester(stk::mesh::MetaData &mesh_meta_data, MPI_Comm comm) :
        stk::mesh::BulkData(mesh_meta_data, comm){}
    virtual ~BulkDataTester() {}

    void my_internal_resolve_shared_modify_delete()
    {
        this->internal_resolve_shared_modify_delete();
    }

    void my_internal_resolve_ghosted_modify_delete()
    {
        this->internal_resolve_ghosted_modify_delete();
    }

    void my_update_comm_list_based_on_changes_in_comm_map()
    {
        this->update_comm_list_based_on_changes_in_comm_map();
    }

    void my_internal_update_distributed_index(std::vector<stk::mesh::Entity> & shared_new )
    {
        this->internal_update_distributed_index( shared_new );
    }

    void my_internal_update_distributed_index(stk::mesh::EntityRank entityRank, std::vector<stk::mesh::Entity> & shared_new )
    {
        this->internal_update_distributed_index(entityRank, shared_new);
    }

    void my_resolve_ownership_of_modified_entities(std::vector<stk::mesh::Entity> & shared_modified )
    {
        this->resolve_ownership_of_modified_entities( shared_modified );
    }

    void my_move_entities_to_proper_part_ownership( std::vector<stk::mesh::Entity> &shared_modified )
    {
        this->move_entities_to_proper_part_ownership( shared_modified );
    }

    void my_update_comm_list( std::vector<stk::mesh::Entity> &shared_modified )
    {
        this->update_comm_list( shared_modified );
    }

    void my_fillLocallyCreatedOrModifiedEntities(stk::parallel::DistributedIndex::KeyTypeVector &local_created_or_modified)
    {
           this->fillLocallyCreatedOrModifiedEntities(local_created_or_modified);
    }

    void my_internal_resolve_shared_membership()
    {
        this->internal_resolve_shared_membership();
    }

    void my_internal_regenerate_aura()
    {
        this->internal_regenerate_aura();
    }

    void reset_closure_count(stk::mesh::Entity entity)
    {
        m_closure_count[entity.local_offset()] = 0;
    }
};

void testElementMove(int fromProc, int toProc, int myProc, int elementToMoveId, stk::mesh::BulkData &stkMeshBulkData);
void populateBulkDataWithFile(const std::string& exodusFileName, MPI_Comm communicator, stk::mesh::BulkData& bulkData);
void checkCommListAndMap(const stk::mesh::BulkData& stkMeshBulkData, bool isAfterIGMD);

void checkStatesOfEntities(std::vector<std::vector<stk::mesh::EntityState> > &nodeStates,
        std::vector<std::vector<stk::mesh::EntityState> > &elementStates,
        bool (&areNodesValid)[2][20], bool (&areElementsValid)[2][4],
        stk::mesh::BulkData &stkMeshBulkData);
void checkThatMeshIsParallelConsistent(stk::mesh::BulkData& stkMeshBulkData);
void mark_element3_as_modified(stk::mesh::BulkData& stkMeshBulkData);
void makeSureEntityIsValidOnCommListAndBulkData(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey);
void makeSureEntityIsValidOnCommListAndBut_NOT_BulkData(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey);
void writeMesh(int myProcId, std::string label, const std::vector<std::string> &meshStart);
void getMeshLineByLine(const stk::mesh::BulkData &stkMeshBulkData, std::vector<std::string> &output);
void fillItemsToFind(std::vector<std::string> &items_to_find);
void getLineNumbersForItemsToFind(const std::vector<std::string>& items_to_find, const std::vector<std::string>& mesh_lines,
        std::vector<int>& line_numbers);
void compareMeshes(int my_proc_id, const std::vector<std::string>& mesh_before, const std::vector<std::string>& mesh_after);
void checkCommMapsAndLists(stk::mesh::BulkData& stkMeshBulkData);
void destroy_element3_on_proc_1(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &elementToDestroyKey);
void checkCommMapsAndListsAfterIRSMD(stk::mesh::BulkData& stkMeshBulkData);
void checkCommMapsAndListsAfterIRGMD(BulkDataTester& stkMeshBulkData);
void create_edges(stk::mesh::BulkData& stkMeshBulkData,
                std::vector<stk::mesh::EntityId>& edgeIds,
                std::vector<std::vector<stk::mesh::EntityId> > &nodeIdsForEdge,
                std::vector<std::vector<stk::mesh::EntityId> > &elementRelations,
                std::vector<stk::mesh::Entity> &edgeEntities, stk::mesh::Part& edge_part);
void checkResultsOfIRSMD_for_edges(stk::mesh::BulkData &stkMeshBulkData);
void checkResultsOfIRGMD_for_edges(BulkDataTester &stkMeshBulkData, std::vector<stk::mesh::Entity> &edgeEntities);
void checkItAllForThisCase(stk::mesh::BulkData &stkMeshBulkData);

TEST(BulkData, ModificationEnd_IR_ghosted_modifiy_delete)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    if(numProcs == 2)
    {
        //============== Load 1x1x4 mesh into BulkData

        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

        std::string exodusFileName = getOption("-i", "generated:1x1x4");
        populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

        //============== Get Global Count Of Elements and Nodes

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);

        size_t numNodes = 20;
        size_t numEdges = 0;
        size_t numFaces = 0;
        size_t numElements = 4;

        ASSERT_EQ(numNodes, globalCounts[stk::topology::NODE_RANK]);
        ASSERT_EQ(numEdges, globalCounts[stk::topology::EDGE_RANK]);
        ASSERT_EQ(numFaces, globalCounts[stk::topology::FACE_RANK]);
        ASSERT_EQ(numElements, globalCounts[stk::topology::ELEMENT_RANK]);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        bool isCheckAfterCallToIRGMD = false;

        checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);
        stkMeshBulkData.modification_begin();
        stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();
        checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);

        std::vector<std::vector<stk::mesh::EntityState> > elementStates(numProcs);
        for (size_t i=0;i<elementStates.size();i++)
        {
            elementStates[i].resize(numElements, stk::mesh::Unchanged);
        }

        std::vector<std::vector<stk::mesh::EntityState> > nodeStates(numProcs);
        for (size_t i=0;i<nodeStates.size();i++)
        {
            nodeStates[i].resize(numNodes, stk::mesh::Unchanged);
        }

        bool areNodesValid[2][20] = {
                { true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, false, false, false, false, },
                { false, false, false, false, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true }
        };

        bool areElementsValid[2][4] = {
                { true, true, true, false, },
                { false, true, true, true }
        };

        checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);

        mark_element3_as_modified(stkMeshBulkData);

        int element3Index = 2;
        int proc_that_modified_element3 = 1;
        elementStates[proc_that_modified_element3][element3Index] = stk::mesh::Modified;

        checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);

        //============== Call IRGMD and check what happens to states of entities and comm maps and lists

        stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();

        int proc_that_no_longer_has_element3 = 0;
        areElementsValid[proc_that_no_longer_has_element3][element3Index] = false;

        for (size_t i=9;i<=16;i++)
        {
            int proc_on_which_nodes_9thru16_are_modified = 0;
            nodeStates[proc_on_which_nodes_9thru16_are_modified][i-1] = stk::mesh::Modified;
        }

        int elementNextToModifiedElementOnProc0 = 1;
        int proc_that_has_element2_modified_as_result_of_element3_modification = 0;
        elementStates[proc_that_has_element2_modified_as_result_of_element3_modification][elementNextToModifiedElementOnProc0] = stk::mesh::Modified;

        checkStatesOfEntities(nodeStates, elementStates, areNodesValid, areElementsValid, stkMeshBulkData);

        isCheckAfterCallToIRGMD = true;
        checkCommListAndMap(stkMeshBulkData, isCheckAfterCallToIRGMD);
    }
}

TEST(BulkData, ModificationEnd)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    if(numProcs == 2)
    {
        //============== Load 1x1x4 mesh into Bulk Data

        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

        // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
        // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12

        std::string exodusFileName = getOption("-i", "generated:1x1x4");
        populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

        //============== Let's destroy element 3 and make a note of node 9 since it's connected between elements 2 and 3

        int elementToDestroy = 3;
        int nodeWhichShouldNoLongerBeSharedOnProc1 = 9;

        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,nodeWhichShouldNoLongerBeSharedOnProc1);
        stk::mesh::Entity node9Entity = stkMeshBulkData.get_entity(nodeEntityKey);
        makeSureEntityIsValidOnCommListAndBulkData(stkMeshBulkData, nodeEntityKey);

        stk::mesh::EntityKey elementToDestroyKey(stk::topology::ELEMENT_RANK, elementToDestroy);
        makeSureEntityIsValidOnCommListAndBulkData(stkMeshBulkData, elementToDestroyKey);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        checkCommMapsAndLists(stkMeshBulkData);

        destroy_element3_on_proc_1(stkMeshBulkData, elementToDestroyKey);

        stkMeshBulkData.my_internal_resolve_shared_modify_delete();

        if ( stkMeshBulkData.parallel_rank() == 0 )
        {
            makeSureEntityIsValidOnCommListAndBulkData(stkMeshBulkData, nodeEntityKey);
        }
        else
        {
            makeSureEntityIsValidOnCommListAndBut_NOT_BulkData(stkMeshBulkData, nodeEntityKey);
        }

        checkCommMapsAndListsAfterIRSMD(stkMeshBulkData);

        //============== Make sure we only have 3 elements remaining

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        EXPECT_EQ(3u, globalCounts[stk::topology::ELEMENT_RANK]);

        stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();

        checkCommMapsAndListsAfterIRGMD(stkMeshBulkData);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        //============== Check the result of internal_update_distributed_index

        std::vector<stk::mesh::Entity> shared_new;
        stkMeshBulkData.my_internal_update_distributed_index( shared_new );

        EXPECT_TRUE(shared_new.empty()==true);
    }
}

TEST(BulkData, ModificationEnd_create_an_edge)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    if(numProcs == 2)
    {
        //============== Load 1x1x4 mesh into Bulk Data

        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

        // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
        // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
        // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

        std::string exodusFileName = getOption("-i", "generated:1x1x4");
        populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

        //============== Before starting, make sure mesh is parallel consistent

        std::vector<std::string> meshStart;
        getMeshLineByLine(stkMeshBulkData, meshStart);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        int myProcId = stkMeshBulkData.parallel_rank();

        std::vector<stk::mesh::EntityId> edgeIds;
        std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
        std::vector<std::vector<stk::mesh::EntityId> > elementRelations;

        edgeIds.push_back(100+myProcId);
        nodeIdsForEdge.resize(edgeIds.size());
        elementRelations.resize(edgeIds.size());
        for ( size_t i=0;i<edgeIds.size();i++)
        {
            nodeIdsForEdge[i].resize(2);
            nodeIdsForEdge[i][0] = 9;
            nodeIdsForEdge[i][1] = 10;

            elementRelations[i].resize(1);
            elementRelations[i][0] = 2+myProcId;
        }

        std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
        stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);

        create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);

        stkMeshBulkData.my_internal_resolve_shared_modify_delete();

        checkResultsOfIRSMD_for_edges(stkMeshBulkData);

        //============== Make sure we have 2 edges

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        size_t numEdgesTotal = 2;
        EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);

        //============== Check results of IRGMD

        stkMeshBulkData.my_internal_resolve_ghosted_modify_delete();

        checkResultsOfIRGMD_for_edges(stkMeshBulkData, edgeEntities);

        //============== Before starting, make sure mesh is parallel consistent

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        //============== Check the result of internal_update_distributed_index

        for (size_t i=0;i<edgeIds.size();i++)
        {
            stk::mesh::EntityKey edgeKey(stk::topology::EDGE_RANK, edgeIds[i]);
            stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), edgeKey);
            bool edge_is_not_in_comm_list = iter == stkMeshBulkData.comm_list().end() || iter->key != edgeKey;
            EXPECT_TRUE(edge_is_not_in_comm_list);

            const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
            EXPECT_FALSE( is_entity_ghosted );

            const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
            EXPECT_FALSE( is_entity_shared );

        }

        for (size_t i=0;i<elementRelations.size();i++)
        {
            for (size_t j=0;j<elementRelations[i].size();j++)
            {
                stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK, elementRelations[i][j]);
                stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
                bool element_not_in_comm_list = iter == stkMeshBulkData.comm_list().end();
                EXPECT_TRUE(element_not_in_comm_list);

                const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
                EXPECT_FALSE( is_entity_ghosted );

                const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
                EXPECT_FALSE( is_entity_shared );
            }
        }

        std::vector<stk::mesh::Entity> shared_new;
        stkMeshBulkData.my_internal_update_distributed_index(stk::topology::EDGE_RANK, shared_new);
        ASSERT_EQ(1u, shared_new.size());
        for (size_t i=0;i<edgeIds.size();i++)
        {
            stk::mesh::EntityId shared_edge_id = 100;
            EXPECT_EQ(shared_edge_id, stkMeshBulkData.identifier(shared_new[i]));

            //============== Funky! Edge is marked as created, not modified. So resolved ownership of "modified" doesn't mess with created.

            EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i]));
            EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(shared_new[i]));
        }

        stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, 9);
        stk::mesh::Entity node = stkMeshBulkData.get_entity(entity_key);

        EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
        }

        stkMeshBulkData.my_resolve_ownership_of_modified_entities(shared_new);

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());

            EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i])) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " failed.";

            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).shared());
            EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));

            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
        }
        EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));

        //============== Hopefully this works!

        stkMeshBulkData.my_move_entities_to_proper_part_ownership( shared_new );

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).shared());
            if ( stkMeshBulkData.parallel_rank() == 0 )
            {
                EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
            }
            else
            {
                EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).owned());
            }
            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));
            EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(shared_new[i])) ;
        }

        stkMeshBulkData.my_update_comm_list( shared_new );

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            stk::mesh::EntityKey edgeKey = stkMeshBulkData.entity_key(shared_new[i]);
            stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), edgeKey);
            bool edge_is_in_comm_list = iter != stkMeshBulkData.comm_list().end() && iter->key == edgeKey;
            EXPECT_TRUE(edge_is_in_comm_list);

            const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
            EXPECT_FALSE( is_entity_ghosted );

            const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
            EXPECT_TRUE( is_entity_shared );
        }

        //============== Not sure what IRSM is supposed to do?

        std::vector<std::string> meshBefore;
        getMeshLineByLine(stkMeshBulkData, meshBefore);

        stkMeshBulkData.my_internal_resolve_shared_membership();

//        std::vector<std::string> meshAfter;
//        getMeshLineByLine(stkMeshBulkData, meshAfter);
//
//        compareMeshes(stkMeshBulkData.parallel_rank(), meshBefore, meshAfter);

        stkMeshBulkData.my_internal_regenerate_aura();

//        std::vector<std::string> meshAfterAura;
//        getMeshLineByLine(stkMeshBulkData, meshAfterAura);
////
////        compareMeshes(stkMeshBulkData.parallel_rank(), meshStart, meshAfterAura);
////        writeMesh(myProcId, "Before::", meshStart);
//        writeMesh(myProcId, "After::", meshAfterAura);
        checkItAllForThisCase(stkMeshBulkData);
    }
}

TEST(BulkData, ModificationEnd_create_an_edge_new)
{
    MPI_Comm communicator = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(communicator, &numProcs);

    if(numProcs == 2)
    {
        //============== Load 1x1x4 mesh into Bulk Data

        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        BulkDataTester stkMeshBulkData(stkMeshMetaData, communicator);

        // Elements 1 and 2 on proc 0, Elements 3 and 4 on proc 1
        // Elements 2 and 3 are shared because of nodes 9, 10, 11, 12
        // Element 2 is ghosted onto Proc 1, and Element 0 is ghosted onto Proc 0

        std::string exodusFileName = getOption("-i", "generated:1x1x4");
        populateBulkDataWithFile(exodusFileName, communicator, stkMeshBulkData);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        int myProcId = stkMeshBulkData.parallel_rank();

        std::vector<stk::mesh::EntityId> edgeIds;
        std::vector<std::vector<stk::mesh::EntityId> > nodeIdsForEdge;
        std::vector<std::vector<stk::mesh::EntityId> > elementRelations;

        edgeIds.push_back(100+myProcId);
        nodeIdsForEdge.resize(edgeIds.size());
        elementRelations.resize(edgeIds.size());
        for ( size_t i=0;i<edgeIds.size();i++)
        {
            nodeIdsForEdge[i].resize(2);
            nodeIdsForEdge[i][0] = 9;
            nodeIdsForEdge[i][1] = 10;

            elementRelations[i].resize(1);
            elementRelations[i][0] = 2+myProcId;
        }

        std::vector<stk::mesh::Entity> edgeEntities(edgeIds.size());
        stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);

        create_edges(stkMeshBulkData, edgeIds, nodeIdsForEdge, elementRelations, edgeEntities, edge_part);

        std::vector<size_t> globalCounts;
        stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
        size_t numEdgesTotal = 2;
        EXPECT_EQ(numEdgesTotal, globalCounts[stk::topology::EDGE_RANK]);

        checkThatMeshIsParallelConsistent(stkMeshBulkData);

        std::vector<stk::mesh::Entity> shared_new;
        stkMeshBulkData.my_internal_update_distributed_index(stk::topology::EDGE_RANK, shared_new);
        ASSERT_EQ(1u, shared_new.size());
        for (size_t i=0;i<edgeIds.size();i++)
        {
            stk::mesh::EntityId shared_edge_id = 100;
            EXPECT_EQ(shared_edge_id, stkMeshBulkData.identifier(shared_new[i]));

            //============== Funky! Edge is marked as created, not modified. So resolved ownership of "modified" doesn't mess with created.

            EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i]));
            EXPECT_EQ(stk::mesh::Created, stkMeshBulkData.state(shared_new[i]));
        }

        stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, 9);
        stk::mesh::Entity node = stkMeshBulkData.get_entity(entity_key);

        EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
        }

        stkMeshBulkData.my_resolve_ownership_of_modified_entities(shared_new);

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());

            EXPECT_EQ(myProcId, stkMeshBulkData.parallel_owner_rank(shared_new[i])) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " failed.";

            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).shared());
            EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));

            EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
        }
        EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));

        //============== Hopefully this works!

        stkMeshBulkData.my_move_entities_to_proper_part_ownership( shared_new );

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).shared());
            if ( stkMeshBulkData.parallel_rank() == 0 )
            {
                EXPECT_TRUE(stkMeshBulkData.bucket(shared_new[i]).owned());
            }
            else
            {
                EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).owned());
            }
            EXPECT_FALSE(stkMeshBulkData.bucket(shared_new[i]).member(stkMeshBulkData.ghosting_part(stkMeshBulkData.aura_ghosting())));
            EXPECT_EQ(0, stkMeshBulkData.parallel_owner_rank(shared_new[i])) ;
        }

        stkMeshBulkData.my_update_comm_list( shared_new );

        for (size_t i=0;i<edgeEntities.size();i++)
        {
            stk::mesh::EntityKey edgeKey = stkMeshBulkData.entity_key(shared_new[i]);
            stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), edgeKey);
            bool edge_is_in_comm_list = iter != stkMeshBulkData.comm_list().end() && iter->key == edgeKey;
            EXPECT_TRUE(edge_is_in_comm_list);

            const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.aura_ghosting()).empty();
            EXPECT_FALSE( is_entity_ghosted );

            const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(edgeKey, stkMeshBulkData.shared_ghosting()).empty();
            EXPECT_TRUE( is_entity_shared );
        }

//        stkMeshBulkData.my_internal_resolve_shared_membership();

//        std::vector<std::string> meshAfterAura;
//        getMeshLineByLine(stkMeshBulkData, meshAfterAura);
//        writeMesh(myProcId, "After::", meshAfterAura);

        checkItAllForThisCase(stkMeshBulkData);
    }
}

void writeMesh(int myProcId, std::string label, const std::vector<std::string> &meshStart)
{
    std::ostringstream msg;
    for (size_t i=0;i<meshStart.size();i++)
    {
        msg.str(std::string());
        msg << "P[" << myProcId << "] " << label << "\t" << meshStart[i] << std::endl;
        std::cerr << msg.str();
    }
}

void getMeshLineByLine(const stk::mesh::BulkData &stkMeshBulkData, std::vector<std::string> &output)
{
    std::ostringstream msg;
    stkMeshBulkData.dump_all_mesh_info(msg);
    std::istringstream iss(msg.str());

    std::string s;
    while ( std::getline(iss, s) )
    {
        output.push_back(s);
    }
}

void fillItemsToFind(std::vector<std::string> &items_to_find)
{
    items_to_find.clear();
    items_to_find.push_back("MetaData info");
    items_to_find.push_back("Entity rank names:");
    items_to_find.push_back("Special Parts:");
    items_to_find.push_back("All parts:");
    items_to_find.push_back("All fields:");
    items_to_find.push_back("BulkData  info");
    items_to_find.push_back("All NODE entities:");
    items_to_find.push_back("All EDGE entities:");
    items_to_find.push_back("All FACE entities:");
    items_to_find.push_back("All ELEMENT entities:");
}

void getLineNumbersForItemsToFind(const std::vector<std::string>& items_to_find, const std::vector<std::string>& mesh_lines,
        std::vector<int>& line_numbers)
{
    int myProcId=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcId);

    line_numbers.clear();
    line_numbers.resize(items_to_find.size()+1);

    size_t foundCounter=0;

    std::ostringstream msg;
    for (size_t i=0;i<mesh_lines.size();i++)
    {
        msg << "[" << myProcId << "] mesh line: " <<  i << "\t" << mesh_lines[i] << std::endl;
        // std::cerr << msg.str();
        msg.str(std::string());
        if ( foundCounter < items_to_find.size() && mesh_lines[i].find(items_to_find[foundCounter]) != std::string::npos )
        {
            line_numbers[foundCounter] = i;
            ++foundCounter;
        }
    }

    ASSERT_EQ(foundCounter, items_to_find.size()) << "could not find " << items_to_find[foundCounter];
    line_numbers[foundCounter]=mesh_lines.size();
}

void compareMeshes(int my_proc_id, const std::vector<std::string>& mesh_before, const std::vector<std::string>& mesh_after)
{
    if ( my_proc_id == 0 )
    {
        std::vector<std::string> items_to_find;
        fillItemsToFind(items_to_find);

        std::vector<int> line_numbers_mesh1;
        getLineNumbersForItemsToFind(items_to_find, mesh_before, line_numbers_mesh1);

        std::vector<int> line_numbers_mesh2;
        getLineNumbersForItemsToFind(items_to_find, mesh_after, line_numbers_mesh2);

        EXPECT_EQ(mesh_before.size(), mesh_after.size());
        ASSERT_EQ(line_numbers_mesh1.size(), line_numbers_mesh2.size());

        for (size_t i=0;i<line_numbers_mesh1.size()-1;i++)
        {
            int numLinesMesh1 = line_numbers_mesh1[i+1]-line_numbers_mesh1[i];
            int numLinesMesh2 = line_numbers_mesh2[i+1]-line_numbers_mesh2[i];
            EXPECT_EQ(numLinesMesh1, numLinesMesh2) <<
                    "[" << my_proc_id << "] Failed for " << items_to_find[i];
            int count = std::min(numLinesMesh1, numLinesMesh2);
            for (int j=0;j<count;j++)
            {
                EXPECT_EQ(mesh_before[line_numbers_mesh1[i]+j],
                          mesh_after[line_numbers_mesh2[i]+j]) << "for proc [" << my_proc_id << "]";
            }
        }
    }
}

void testElementMove(int fromProc, int toProc, int myProc, int elementToMoveId, stk::mesh::BulkData &stkMeshBulkData)
{
    stk::mesh::Entity elementToMove = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToMoveId);
    std::vector<std::pair<stk::mesh::Entity, int> > entityProcPairs;
    if(myProc == fromProc)
    {
        entityProcPairs.push_back(std::make_pair(elementToMove, toProc));
    }
    stkMeshBulkData.change_entity_owner(entityProcPairs);
}

void populateBulkDataWithFile(const std::string& exodusFileName, MPI_Comm communicator, stk::mesh::BulkData& bulkData)
// STK IO module will be described in separate chapter.
// It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
// The order of the following lines in {} are important
{
  stk::io::StkMeshIoBroker exodusFileReader(communicator);

  // Inform STK IO which STK Mesh objects to populate later
  exodusFileReader.set_bulk_data(bulkData);

  exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

  // Populate the MetaData which has the descriptions of the Parts and Fields.
  exodusFileReader.create_input_mesh();

  // Populate entities in STK Mesh from Exodus file
  exodusFileReader.populate_bulk_data();
}

void checkCommListAndMap(const stk::mesh::BulkData& stkMeshBulkData, bool isAfterIRGMD)
{
    bool doesElementAppearInAuraCommMap[2][4] = { { false, true, true, false, },
                                               { false, true, true, false } };

    bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                            { false, true, true, false } };

    if ( isAfterIRGMD )
    {
        doesElementAppearInAuraCommMap[0][1] = false;
        doesElementAppearInAuraCommMap[0][2] = false;
        doesElementAppearInAuraCommMap[1][2] = false;
        isElementValidInCommList[0][2] = false;
        isElementValidInCommList[0][2] = false;
    }

    stk::mesh::EntityCommListInfoVector::const_iterator iter;
    int myProcId = stkMeshBulkData.parallel_rank();
    for (unsigned int element=1;element<=4;element++)
    {
        int elementoffset = element-1;
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        if ( iter == stkMeshBulkData.comm_list().end() || iter->key != elementKey )
        {
            EXPECT_EQ(isElementValidInCommList[myProcId][elementoffset], false) <<
                    "Proc " << myProcId << " for element " << element;
        }
        else
        {
            EXPECT_EQ(isElementValidInCommList[myProcId][elementoffset], stkMeshBulkData.is_valid(iter->entity)) <<
                    "proc " << myProcId << " using offset " << elementoffset;
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ(doesElementAppearInAuraCommMap[myProcId][elementoffset], is_entity_ghosted ) <<
                "proc " << myProcId << " using offset " << elementoffset;
    }
}


void checkStatesOfEntities(std::vector<std::vector<stk::mesh::EntityState> > &nodeStates,
        std::vector<std::vector<stk::mesh::EntityState> > &elementStates,
        bool (&areNodesValid)[2][20], bool (&areElementsValid)[2][4],
        stk::mesh::BulkData &stkMeshBulkData)
{
    int myProcId = stkMeshBulkData.parallel_rank();
    for (unsigned nodeId=1;nodeId<=nodeStates[myProcId].size();nodeId++)
    {
        stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, nodeId);
        stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
        ASSERT_EQ(areNodesValid[myProcId][nodeId-1], stkMeshBulkData.is_valid(entity)) <<
                "Proc " << myProcId << " using node " << nodeId;
        if ( stkMeshBulkData.is_valid(entity) )
        {
            EXPECT_EQ(nodeStates[myProcId][nodeId-1], stkMeshBulkData.state(entity));
        }
    }

    for (unsigned elementId=1;elementId<=elementStates[myProcId].size();elementId++)
    {
        stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, elementId);
        stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
        ASSERT_EQ(areElementsValid[myProcId][elementId-1], stkMeshBulkData.is_valid(entity));
        if ( stkMeshBulkData.is_valid(entity) )
        {
            EXPECT_EQ(elementStates[myProcId][elementId-1], stkMeshBulkData.state(entity)) ;
        }
    }
}

void mark_element3_as_modified(stk::mesh::BulkData& stkMeshBulkData)
{
    int elementToModify = 3;
    stk::mesh::EntityKey elementToModifyKey(stk::topology::ELEMENT_RANK, elementToModify);
    stk::mesh::Entity entity = stkMeshBulkData.get_entity(elementToModifyKey);

    if ( stkMeshBulkData.parallel_rank() == 1 )
    {
        stkMeshBulkData.set_state(entity, stk::mesh::Modified);
    }
}

stk::mesh::EntityCommListInfoVector::const_iterator makeSureEntityIsValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), entityKey);
    EXPECT_TRUE(iter != stkMeshBulkData.comm_list().end());
    EXPECT_EQ(entityKey,iter->key);
    return iter;
}

stk::mesh::EntityCommListInfoVector::const_iterator makeSureEntityIs_NOT_ValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), entityKey);
    EXPECT_TRUE(iter == stkMeshBulkData.comm_list().end() || entityKey != iter->key);
    return iter;
}

void makeSureEntityIsValidOnCommListAndBut_NOT_BulkData(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
    stk::mesh::EntityCommListInfoVector::const_iterator iter = makeSureEntityIsValidOnCommList(stkMeshBulkData, entityKey);
    EXPECT_FALSE(stkMeshBulkData.is_valid(iter->entity));
}

void makeSureEntityIsValidOnCommListAndBulkData(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
    stk::mesh::EntityCommListInfoVector::const_iterator iter = makeSureEntityIsValidOnCommList(stkMeshBulkData, entityKey);
    EXPECT_TRUE(stkMeshBulkData.is_valid(iter->entity));
}

void checkThatMeshIsParallelConsistent(stk::mesh::BulkData& stkMeshBulkData)
{
    std::ostringstream msg ;
    bool is_consistent = true;
    is_consistent = comm_mesh_verify_parallel_consistency( stkMeshBulkData , msg );
    EXPECT_TRUE(is_consistent) << msg.str();
}

void checkCommMapsAndLists(stk::mesh::BulkData& stkMeshBulkData)
{
    bool isNodeInCommList[2][20] = {
            {false, false, false, false,
             true,  true,  true,  true,
             true,  true,  true,  true,
             true,  true,  true,  true,
             false, false, false, false,},

            {false, false, false, false,
             true,  true,  true,  true,
             true,  true,  true,  true,
             true,  true,  true,  true,
             false, false, false, false}
    };

    bool isNodeInGhostedCommMap[2][20] = {
            {false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
            },
            {false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
            }
    };

    bool IsNodeInSharedCommMap[2][20] = {
            {false, false, false, false,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             false, false, false, false,
            },
            {false, false, false, false,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             false, false, false, false,
            }
    };

    for (unsigned int nodeId=1;nodeId<=20;nodeId++)
    {
        stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), nodeKey);

        if ( iter != stkMeshBulkData.comm_list().end() && iter->key == nodeKey )
        {
            EXPECT_TRUE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
        }
        else
        {
            EXPECT_FALSE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ( isNodeInGhostedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_ghosted ) <<
                "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_EQ( IsNodeInSharedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_shared ) <<
                "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }
}

void destroy_element3_on_proc_1(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::EntityKey &elementToDestroyKey)
{
    stkMeshBulkData.modification_begin();

    ASSERT_TRUE ( stkMeshBulkData.is_valid(stkMeshBulkData.get_entity(elementToDestroyKey)) );

    if ( stkMeshBulkData.parallel_rank() == 1 )
    {
        stkMeshBulkData.destroy_entity( stkMeshBulkData.get_entity(elementToDestroyKey) );
    }
}

void checkCommMapsAndListsAfterIRSMD(stk::mesh::BulkData& stkMeshBulkData)
{
    //============== Checking result of IRSMD


    /*
     * We deleted element 3 from mesh that looks like E1 : E2 : E3 : E4
     * Then we called internal_resolve_shared_modify_delete. Note,
     * nodes 9 thru 12 are on the boundary between E2 and E3. Also, proc 0
     * has E1 : E2 and proc 1 has E3 : E4. So the proc decomp looks like
     *
     * proc 0                               proc 1
     *    E1 : E2 : GE3                             GE2 : E3 : E4
     *
     * where GE2 and GE3 are ghosted elements.
     *
     * Result:
     *      For proc 0: nodes 9 thru 12 are still valid entities
     *      For proc 1: nodes 9 thru 12 are not valid anymore
     *
     *      Proc 0 still *thinks* it needes to communicate about elements 1, 2, and 3. Why?
     *      Proc 1 doesn't have any elements it needs to communicate with.
     *
     *      Proc 0 believes elements 2 and 3 are still needed for communication for aura.
     *      Proc 1 believes elements 2 and 3 are still needed for communication for aura.
     *
     * Question: Ghosting is NOT right, but is sharing right?
     *      According to test below, nothing is shared. So yes!?
     */

    bool isNodeValidInCommList[2][4] = { { true, true, true, true, },
                                         { false, false, false, false } };

    bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                            { false, false, false, false } };

    bool isElementInAuraCommMap[2][4] = { { false, true, true, false, },
                                             { false, true, true, false } };

    //============== Testing against the boolean vectors

    for (unsigned int node=9;node<=12;node++)
    {
        int nodeoffset = node-9;
        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,node);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), nodeEntityKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != nodeEntityKey)
        {
            EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], false);
        }
        else
        {
            EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(nodeEntityKey, stkMeshBulkData.aura_ghosting()).empty();
        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(nodeEntityKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
        EXPECT_FALSE( is_entity_ghosted );
    }

    for (unsigned int element=1;element<=4;element++)
    {
        int elementoffset = element-1;
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != elementKey)
        {
            EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], false);
        }
        else
        {
            EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ(isElementInAuraCommMap[stkMeshBulkData.parallel_rank()][elementoffset], is_entity_ghosted);

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }
}

void checkCommMapsAndListsAfterIRGMD(BulkDataTester& stkMeshBulkData)
{
    bool isElementValidInCommListAfterIRGMD[2][4] = { { false, true, false, false, },
                                                      { false, false, false, false } };

    //============== Check results against boolean above

    for (unsigned int element=1;element<=4;element++)
    {
        int elementoffset = element-1;
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != elementKey)
        {
            EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], false);
        }
        else
        {
            EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_FALSE( is_entity_ghosted );

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }

    //============== Check results of updating the comm list based on changes in the comm map

    stkMeshBulkData.my_update_comm_list_based_on_changes_in_comm_map();

    //============== No element is ghosted, No node is shared, and comm list is empty

    for (unsigned int element=1;element<=4;element++)
    {
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        EXPECT_TRUE(iter == stkMeshBulkData.comm_list().end());

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_FALSE( is_entity_ghosted );

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }

    // Element 3 has been delete so:
    for (unsigned int nodeId=1;nodeId<=20;nodeId++)
    {
        stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), nodeKey);

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
        if ( nodeId >=5 and nodeId <=8)
        {
            EXPECT_TRUE( iter != stkMeshBulkData.comm_list().end() && iter->key == nodeKey );
            EXPECT_TRUE( is_entity_ghosted );
        }
        else
        {
            EXPECT_FALSE ( iter != stkMeshBulkData.comm_list().end() && iter->key == nodeKey ) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

            EXPECT_FALSE( is_entity_ghosted ) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
        }

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared ) <<
                "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }
}

void connectElementToEdge(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity element,
        stk::mesh::Entity edge, stk::mesh::Permutation perm, stk::mesh::OrdinalVector& ordinal_scratch,
        stk::mesh::PartVector& part_scratch, const std::vector<stk::mesh::EntityId> nodeIdsForEdge)
{
    ASSERT_TRUE(stkMeshBulkData.is_valid(element));
    stk::topology elem_top = stkMeshBulkData.bucket(element).topology();
    std::vector<stk::mesh::EntityId> nodeIds(2);
    std::vector<stk::mesh::EntityId> elementNodes(elem_top.num_nodes());
    unsigned edge_ordinal = 100000;
    stk::mesh::Entity const * elem_nodes = stkMeshBulkData.begin_nodes(element);
    for (size_t k=0;k<elem_top.num_nodes();k++)
    {
        elementNodes[k] =stkMeshBulkData.identifier(elem_nodes[k]);
    }

    for (size_t k=0;k<elem_top.num_edges();k++)
    {
        elem_top.edge_nodes(elementNodes, k, nodeIds.begin());
        stk::topology edge_top = elem_top.edge_topology();
        if ( edge_top.equivalent(nodeIds, nodeIdsForEdge).first )
        {
            edge_ordinal = k;
            break;
        }
    }
    ASSERT_TRUE(edge_ordinal != 100000);
    stkMeshBulkData.declare_relation(element, edge, edge_ordinal, perm, ordinal_scratch, part_scratch);
}

void create_edges(stk::mesh::BulkData& stkMeshBulkData, std::vector<stk::mesh::EntityId>& edgeIds,
                std::vector<std::vector<stk::mesh::EntityId> > &nodeIdsForEdge,
                std::vector<std::vector<stk::mesh::EntityId> > &elementRelations,
                std::vector<stk::mesh::Entity> &edgeEntities,
                stk::mesh::Part& edge_part)
{
    //============== Create one edge between elements 2 and 3 (nodes 9 and 10) with edge id 100

    stk::mesh::OrdinalVector ordinal_scratch;
    ordinal_scratch.reserve(64);
    stk::mesh::PartVector part_scratch;
    part_scratch.reserve(64);

    stk::mesh::PartVector add_parts;
    add_parts.push_back( & stkMeshBulkData.mesh_meta_data().get_cell_topology_root_part( stk::mesh::get_cell_topology( stk::topology::LINE_2 )));
    add_parts.push_back(&edge_part);
    stkMeshBulkData.modification_begin();

    for (size_t i=0;i<edgeIds.size();i++)
    {
        stk::mesh::Entity edge = stkMeshBulkData.declare_entity( stk::topology::EDGE_RANK, edgeIds[i], add_parts);
        edgeEntities[i] = edge;

        std::vector<stk::mesh::Entity> ghostedElements(10);
        stk::mesh::Permutation perm = static_cast<stk::mesh::Permutation>(0);
        {
            std::vector<stk::mesh::Entity> nodes(2);
            ASSERT_EQ(2u, nodeIdsForEdge[i].size());
            for (size_t n=0; n<nodeIdsForEdge[i].size(); ++n)
            {
                stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,nodeIdsForEdge[i][n]);
                stk::mesh::Entity node = stkMeshBulkData.get_entity(nodeEntityKey);
                ASSERT_TRUE(stkMeshBulkData.is_valid(node));
                stkMeshBulkData.declare_relation(edge, node,n, perm, ordinal_scratch, part_scratch);
                EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
                nodes[n] = node;
            }
            stk::mesh::Entity const * elemStartNode1 = stkMeshBulkData.begin_elements(nodes[0]);
            stk::mesh::Entity const * elemEndNode1 = stkMeshBulkData.end_elements(nodes[0]);
            stk::mesh::Entity const * elemStartNode2 = stkMeshBulkData.begin_elements(nodes[1]);
            stk::mesh::Entity const * elemEndNode2 = stkMeshBulkData.end_elements(nodes[1]);

            std::vector<stk::mesh::Entity> elems1(elemStartNode1, elemEndNode1);
            std::sort(elems1.begin(), elems1.end());
            std::vector<stk::mesh::Entity> elems2(elemStartNode2, elemEndNode2);
            std::sort(elems2.begin(), elems2.end());

            std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( elems1.begin(), elems1.end(),
                    elems2.begin(), elems2.end(), ghostedElements.begin());

            ghostedElements.resize(iter-ghostedElements.begin());
        }

        for (size_t j=0;j<ghostedElements.size();j++)
        {
            if ( !stkMeshBulkData.bucket(ghostedElements[j]).owned() )
            {
                connectElementToEdge(stkMeshBulkData, ghostedElements[j], edge, perm, ordinal_scratch, part_scratch, nodeIdsForEdge[i]);
            }
        }

        for (size_t j=0;j<elementRelations[i].size();j++)
        {
            stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,elementRelations[i][j]);
            stk::mesh::Entity element = stkMeshBulkData.get_entity(elementKey);
            connectElementToEdge(stkMeshBulkData, element, edge, perm, ordinal_scratch, part_scratch, nodeIdsForEdge[i]);
        }
    }
}

void checkResultsOfIRSMD_for_edges(stk::mesh::BulkData &stkMeshBulkData)
{
    bool isNodeValidInCommList[2][4] = { { true, true, true, true, },
                                         { true, true, true, true } };

    //  bool isEdgeInCommList[2][1] = { {true}, {false} };
    //

    bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                            { false, true, true, false } };

    bool isElementInAuraCommMap[2][4] = { { false, true, true, false, },
                                          { false, true, true, false } };

    stk::mesh::EntityCommListInfoVector::const_iterator iter;

    //============== Testing against the boolean vectors

    for (unsigned int node=9;node<=12;node++)
    {
        int nodeoffset = node-9;
        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,node);
        iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), nodeEntityKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != nodeEntityKey)
        {
            EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], false);
        }
        else
        {
            EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(nodeEntityKey, stkMeshBulkData.aura_ghosting()).empty();
        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(nodeEntityKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_TRUE( is_entity_shared );
        EXPECT_FALSE( is_entity_ghosted );
    }

    for (unsigned int element=1;element<=4;element++)
    {
        int elementoffset = element-1;
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != elementKey)
        {
            EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], false);
        }
        else
        {
            EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ(isElementInAuraCommMap[stkMeshBulkData.parallel_rank()][elementoffset], is_entity_ghosted);

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }
}

void checkResultsOfIRGMD_for_edges(BulkDataTester &stkMeshBulkData, std::vector<stk::mesh::Entity> &edgeEntities)
{
    for (size_t i=0;i<edgeEntities.size();i++)
    {
        EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
    }

    bool isElementValidInCommListAfterIRGMD[2][4] = { { false, true, false, false, },
                                                      { false, false, true, false } };

    //============== Check results against boolean above

    for (unsigned int element=1;element<=4;element++)
    {
        int elementoffset = element-1;
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        if (iter == stkMeshBulkData.comm_list().end() || iter->key != elementKey)
        {
            EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], false);
        }
        else
        {
            EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_FALSE( is_entity_ghosted );

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }

    //============== Check results of updating the comm list based on changes in the comm map

    stkMeshBulkData.my_update_comm_list_based_on_changes_in_comm_map();

    //============== No element is ghosted, No node is shared, and comm list is empty

    for (unsigned int element=1;element<=4;element++)
    {
        stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), elementKey);
        EXPECT_TRUE(iter == stkMeshBulkData.comm_list().end());

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_FALSE( is_entity_ghosted );

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_FALSE( is_entity_shared );
    }

    bool isNodeInCommList[2][20] = {
            {false, false, false, false,
             true,  true,  true,  true,
             true,  true,  true,  true,
             true,  true,  true,  true,
             false, false, false, false,},

            {false, false, false, false,
             true,  true,  true,  true,
             true,  true,  true,  true,
             true,  true,  true,  true,
             false, false, false, false}
    };

    bool isNodeInGhostedCommMap[2][20] = {
            {false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
            },
            {false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
            }
    };

    bool IsNodeInSharedCommMap[2][20] = {
            {false, false, false, false,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             false, false, false, false,
            },
            {false, false, false, false,
             false, false, false, false,
             true,  true,  true,  true,
             false, false, false, false,
             false, false, false, false,
            }
    };

    for (unsigned int nodeId=1;nodeId<=20;nodeId++)
    {
        stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
        stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), nodeKey);

        if ( iter != stkMeshBulkData.comm_list().end() && iter->key == nodeKey )
        {
            EXPECT_TRUE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
        }
        else
        {
            EXPECT_FALSE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                    "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
        }

        const bool is_entity_ghosted = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ( isNodeInGhostedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_ghosted ) <<
                "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

        const bool is_entity_shared = !stkMeshBulkData.entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_EQ( IsNodeInSharedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_shared ) <<
                "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }
}

void checkResults(stk::mesh::BulkData& stkMeshBulkData,
        const size_t numEntities,
        const stk::mesh::Part& edge_part,
        const std::vector<stk::mesh::EntityId>& entityIds,
        stk::mesh::EntityRank rank,
        bool *isEntityValidOnBulkData,
        bool *isEntityValidOnCommList,
        bool *isEntityOwned,
        bool *isEntityGhosted,
        bool *isEntityShared,
        bool *isEntityOnEdgePart,
        bool *isEntityOnAuraCommMap,
        bool *isEntityOnSharedCommMap)
{
    int procId=-1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);

    for (size_t i=0;i<numEntities;i++)
    {
        stk::mesh::EntityKey entity_key(rank, entityIds[i]);
        stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
        EXPECT_EQ(isEntityValidOnBulkData[i], stkMeshBulkData.is_valid(entity)) << "P[" << procId << "] for rank: " << rank;
        if ( isEntityValidOnCommList[i] )
        {
            makeSureEntityIsValidOnCommList(stkMeshBulkData, entity_key);
        }
        else
        {
            makeSureEntityIs_NOT_ValidOnCommList(stkMeshBulkData, entity_key);
        }
        if ( isEntityValidOnBulkData[i] )
        {
            EXPECT_EQ(isEntityOwned[i], stkMeshBulkData.bucket(entity).owned())<< "P[" << procId << "] for rank: " << rank;
            EXPECT_EQ(isEntityGhosted[i], stkMeshBulkData.bucket(entity).member(stkMeshBulkData.mesh_meta_data().aura_part()))<< "P[" << procId << "] for rank: " << rank;
            EXPECT_EQ(isEntityShared[i], stkMeshBulkData.bucket(entity).shared())<< "P[" << procId << "] for rank: " << rank;
            EXPECT_EQ(isEntityOnEdgePart[i], stkMeshBulkData.bucket(entity).member(edge_part))<< "P[" << procId << "] for rank: " << rank;
        }
        else
        {
            EXPECT_FALSE(isEntityOwned[i])<< "P[" << procId << "] for rank: " << rank;
            EXPECT_FALSE(isEntityGhosted[i])<< "P[" << procId << "] for rank: " << rank;
            EXPECT_FALSE(isEntityShared[i])<< "P[" << procId << "] for rank: " << rank;
            EXPECT_FALSE(isEntityOnEdgePart[i])<< "P[" << procId << "] for rank: " << rank;
        }
        const bool is_node_on_aura_comm_map = !stkMeshBulkData.entity_comm_map(entity_key, stkMeshBulkData.aura_ghosting()).empty();
        EXPECT_EQ( isEntityOnAuraCommMap[i], is_node_on_aura_comm_map )<< "P[" << procId << "] for rank: " << rank;

        const bool is_node_on_shared_comm_map = !stkMeshBulkData.entity_comm_map(entity_key, stkMeshBulkData.shared_ghosting()).empty();
        EXPECT_EQ( isEntityOnSharedCommMap[i], is_node_on_shared_comm_map )<< "P[" << procId << "] for rank: " << rank;
    }
}

void checkEntityRelations(int procId, stk::mesh::BulkData& stkMeshBulkData)
{
    int counter=0;
    {
        stk::mesh::EntityId elementConn[24] = {
                1, 2, 4, 3, 5, 6, 8, 7,
                5, 6, 8, 7, 9, 10, 12, 11,
                9, 10, 12, 11, 13, 14, 16, 15
        };

        for (size_t i=0;i<3;i++)
        {
            stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, i+1+procId);
            stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
            stk::mesh::Entity const * nodesbegin = stkMeshBulkData.begin_nodes(entity);
            ASSERT_TRUE(nodesbegin!=0) << "for proc " << procId;
            stk::mesh::Entity const * nodesend = stkMeshBulkData.end_nodes(entity);
            for (stk::mesh::Entity const * node = nodesbegin; node != nodesend; ++node)
            {
                stk::mesh::EntityKey node_key(stk::topology::NODE_RANK, elementConn[counter]+4*procId);
                stk::mesh::Entity conn_node = stkMeshBulkData.get_entity(node_key);
                EXPECT_EQ(*node, conn_node) << "for proc " << procId;
                counter++;
            }
        }
    }

    stk::mesh::EntityKey edge_key(stk::topology::EDGE_RANK, 100);
    stk::mesh::Entity edge = stkMeshBulkData.get_entity(edge_key);

    stk::mesh::EntityKey element2_key(stk::topology::ELEMENT_RANK, 2);
    stk::mesh::Entity element2 = stkMeshBulkData.get_entity(element2_key);

    stk::mesh::EntityKey element3_key(stk::topology::ELEMENT_RANK, 3);
    stk::mesh::Entity element3 = stkMeshBulkData.get_entity(element3_key);

    {
        for (size_t i=2;i<=3;i++)
        {
            stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, i);
            stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
            stk::mesh::Entity const * edges_begin = stkMeshBulkData.begin_edges(entity);
            ASSERT_TRUE(edges_begin!=0) << "for proc " << procId << " against element " << entity_key;
            EXPECT_EQ( *edges_begin, edge) << "for proc " << procId << " against element " << entity_key;;
        }
    }

    stk::mesh::EntityKey node9_key(stk::topology::NODE_RANK, 9);
    stk::mesh::Entity node9 = stkMeshBulkData.get_entity(node9_key);

    stk::mesh::EntityKey node10_key(stk::topology::NODE_RANK, 10);
    stk::mesh::Entity node10 = stkMeshBulkData.get_entity(node10_key);

    {
        stk::mesh::Entity const * entity = 0;
        entity = stkMeshBulkData.begin_nodes(edge);
        ASSERT_TRUE(entity!=0) << "for proc " << procId;

        EXPECT_EQ(*entity, node9) << "for proc " << procId;
        entity++;
        EXPECT_EQ(*entity, node10) << "for proc " << procId;
    }

    {
        stk::mesh::Entity const * elements_begin = stkMeshBulkData.begin_elements(edge);
        ASSERT_TRUE(elements_begin!=0) << "for proc " << procId;

        bool one_of_elements = (*elements_begin == element2) || (*elements_begin == element3);
        EXPECT_TRUE(one_of_elements) << "for proc " << procId;
        elements_begin++;
        one_of_elements = *elements_begin == element2 || *elements_begin == element3;
        EXPECT_TRUE(one_of_elements) << "for proc " << procId;
    }

    {
        stk::mesh::Entity const * entity = 0;

        entity = stkMeshBulkData.begin_edges(node9);
        ASSERT_TRUE(entity!=0) << "for proc " << procId;
        EXPECT_EQ(*entity, edge) << "for proc " << procId;

        entity = stkMeshBulkData.begin_edges(node10);
        ASSERT_TRUE(entity!=0) << "for proc " << procId;
        EXPECT_EQ(*entity, edge) << "for proc " << procId;

        entity = stkMeshBulkData.begin_elements(node9);
        ASSERT_TRUE(entity!=0) << "for proc " << procId;
        EXPECT_EQ(*entity, element3) << "for proc " << procId;
        entity++;
        EXPECT_EQ(*entity, element2) << "for proc " << procId;

        entity = stkMeshBulkData.begin_elements(node10);
        ASSERT_TRUE(entity!=0) << "for proc " << procId;
        EXPECT_EQ(*entity, element3) << "for proc " << procId;
        entity++;
        EXPECT_EQ(*entity, element2) << "for proc " << procId;
    }
}
void check_it_all_for_proc_0(stk::mesh::BulkData &stkMeshBulkData)
{
    size_t numNodes = 20;
    bool isNodeValidOnBulkData[20] = {
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            false, false, false, false
    };
    bool isNodeValidOnCommList[20] = {
            false, false, false, false,
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            false, false, false, false
    };
    bool isNodeOwned[20] = {
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            false, false, false, false,
            false, false, false, false
    };
    bool isNodeShared[20] = {
            false, false, false, false,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            false, false, false, false
    };
    bool isNodeGhosted[20] = {
            false, false, false, false,
            false, false, false, false,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false
    };
    bool isNodeOnAuraCommMap[20] = {
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false
    };
    bool isNodeOnSharedCommMap[20] = {
            false, false, false, false,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            false, false, false, false
    };
    bool isNodeOnEdgePart[20] = {
                false, false, false, false,
                false, false, false, false,
                true, true, false, false,
                false, false, false, false,
                false, false, false, false
    };
    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
    std::vector<stk::mesh::EntityId> nodeIds(numNodes);
    for (size_t i=0;i<nodeIds.size();i++)
    {
        nodeIds[i] = i+1;
    }
    checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
                 isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

    size_t numEdges = 1;
    bool isEdgeValidOnBulkData[1] = {
            true
    };
    bool isEdgeValidOnCommList[1] = {
            true
    };
    bool isEdgeOwned[1] = {
            true
    };
    bool isEdgeShared[1] = {
                true
    };
    bool isEdgeGhosted[1] = {
            false
    };
    bool isEdgeOnAuraCommMap[1] = {
            false
    };
    bool isEdgeOnSharedCommMap[1] = {
            true
    };
    bool isEdgeOnEdgePart[1] = {
            true
    };
    std::vector<stk::mesh::EntityId> edgeIds(numEdges);
    edgeIds[0] = 100;
    checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
                 isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

    size_t numElements = 4;
    bool isElementValidOnBulkData[4] = {
            true, true, true, false
    };
    bool isElementValidOnCommList[4] = {
            false, true, true, false
    };
    bool isElementOwned[4] = {
            true, true, false, false
    };
    bool isElementShared[4] = {
            false, false, false, false
    };
    bool isElementGhosted[4] = {
            false, false, true, false
    };
    bool isElementOnAuraCommMap[4] = {
            false, true, true, false
    };
    bool isElementOnSharedCommMap[4] = {
            false, false, false, false
    };
    bool isElementOnEdgePart[4] = {
            false, false, false, false
    };
    std::vector<stk::mesh::EntityId> elementIds(numElements);
    for (size_t i=0;i<elementIds.size();i++)
    {
        elementIds[i] = i+1;
    }
    checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
                 isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

    checkEntityRelations(0, stkMeshBulkData);
}

void check_it_all_for_proc_1(stk::mesh::BulkData &stkMeshBulkData)
{
    size_t numNodes = 20;
    bool isNodeValidOnBulkData[20] = {
            false, false, false, false,
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            true, true, true, true
    };

    bool isNodeValidOnCommList[20] = {
            false, false, false, false,
            true, true, true, true,
            true, true, true, true,
            true, true, true, true,
            false, false, false, false
    };

    bool isNodeOwned[20] = {
            false, false, false, false,
            false, false, false, false,
            false, false, false, false,
            true, true, true, true,
            true, true, true, true
    };
    bool isNodeShared[20] = {
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false
    };
    bool isNodeGhosted[20] = {
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            false, false, false, false,
            false, false, false, false
    };
    bool isNodeOnAuraCommMap[20] = {
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false
    };
    bool isNodeOnSharedCommMap[20] = {
            false, false, false, false,
            false, false, false, false,
            true, true, true, true,
            false, false, false, false,
            false, false, false, false
    };
    bool isNodeOnEdgePart[20] = {
            false, false, false, false,
            false, false, false, false,
            true, true, false, false,
            false, false, false, false,
            false, false, false, false
    };
    stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
    std::vector<stk::mesh::EntityId> nodeIds(numNodes);
    for (size_t i=0;i<nodeIds.size();i++)
    {
        nodeIds[i] = i+1;
    }
    checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
                 isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

    size_t numEdges = 1;
    bool isEdgeValidOnBulkData[1] = {
            true
    };
    bool isEdgeValidOnCommList[1] = {
            true
    };
    bool isEdgeOwned[1] = {
            false
    };
    bool isEdgeShared[1] = {
                true
    };
    bool isEdgeGhosted[1] = {
            false
    };
    bool isEdgeOnAuraCommMap[1] = {
            false
    };
    bool isEdgeOnSharedCommMap[1] = {
            true
    };
    bool isEdgeOnEdgePart[1] = {
            true
    };
    std::vector<stk::mesh::EntityId> edgeIds(numEdges);
    edgeIds[0] = 100;
    checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
                 isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

    size_t numElements = 4;
    bool isElementValidOnBulkData[4] = {
            false, true, true, true
    };
    bool isElementValidOnCommList[4] = {
            false, true, true, false
    };
    bool isElementOwned[4] = {
            false, false, true, true
    };
    bool isElementShared[4] = {
            false, false, false, false
    };
    bool isElementGhosted[4] = {
            false, true, false, false
    };
    bool isElementOnAuraCommMap[4] = {
            false, true, true, false
    };
    bool isElementOnSharedCommMap[4] = {
            false, false, false, false
    };
    bool isElementOnEdgePart[4] = {
            false, false, false, false
    };
    std::vector<stk::mesh::EntityId> elementIds(numElements);
    for (size_t i=0;i<elementIds.size();i++)
    {
        elementIds[i] = i+1;
    }
    checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
                 isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

    checkEntityRelations(1, stkMeshBulkData);
}

void checkItAllForThisCase(stk::mesh::BulkData &stkMeshBulkData)
{
    std::vector<size_t> globalCounts;
    stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
    EXPECT_EQ(20u, globalCounts[stk::topology::NODE_RANK]);
    EXPECT_EQ(1u, globalCounts[stk::topology::EDGE_RANK]);
    EXPECT_EQ(0u, globalCounts[stk::topology::FACE_RANK]);
    EXPECT_EQ(4u, globalCounts[stk::topology::ELEMENT_RANK]);

    if ( stkMeshBulkData.parallel_rank() == 0)
    {
        check_it_all_for_proc_0(stkMeshBulkData);
    }
    else
    {
        check_it_all_for_proc_1(stkMeshBulkData);
    }
}

} // end namespace
