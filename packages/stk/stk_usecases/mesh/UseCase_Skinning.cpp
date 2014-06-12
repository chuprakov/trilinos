/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <mesh/UseCase_Skinning.hpp>

#include <stk_mesh/base/BulkModification.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityCommDatabase.hpp>

#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/BoundaryAnalysis.hpp>
#include <stk_mesh/base/SkinMesh.hpp>

#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace {

static const stk::mesh::EntityRank NODE_RANK = stk::topology::NODE_RANK;

void find_owned_nodes_with_relations_outside_closure(
    const stk::mesh::BulkData &mesh,
    stk::mesh::EntityVector & closure,
    stk::mesh::Selector       select_owned,
    stk::mesh::EntityVector & nodes)
{
  nodes.clear();

  stk::mesh::EntityLess lesser(mesh);

  //the closure is a sorted unique vector
  const stk::mesh::EntityRank upward_rank = stk::topology::EDGE_RANK;
  const stk::mesh::EntityId base_id = 0;
  stk::mesh::EntityVector::iterator node_end = std::lower_bound(closure.begin(),
      closure.end(),
      stk::mesh::EntityKey(upward_rank, base_id),
      lesser);

  for (stk::mesh::EntityVector::iterator itr = closure.begin(); itr != node_end; ++itr) {
    stk::mesh::Entity node = *itr;

    if (select_owned(mesh.bucket(node))) {

      bool pushed_back = false;

      for (stk::mesh::EntityRank irank = stk::topology::BEGIN_RANK;
             !pushed_back && (irank != stk::topology::END_RANK);
             ++irank)
      {
        //loop over the relations and check to see if they are in the closure
        stk::mesh::Entity const * relations_iter = mesh.begin(node, irank);
        stk::mesh::Entity const * relations_end = mesh.end(node, irank);

        for (; relations_iter != relations_end; ++relations_iter) {
          stk::mesh::Entity current_entity = *relations_iter;

          //has relation outside of closure
          if ( !std::binary_search(node_end,
                                   closure.end(),
                                   current_entity,
                                   lesser) )
          {
            nodes.push_back(node);
            pushed_back = false;
            break;
          }
        }
      }
    }
  }
}

void copy_nodes_and_break_relations( stk::mesh::BulkData     & mesh,
                                     stk::mesh::EntityVector & closure,
                                     stk::mesh::EntityVector & nodes,
                                     stk::mesh::EntityVector & new_nodes)
{
  stk::mesh::EntityLess lesser(mesh);

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk::mesh::Entity entity = nodes[i];
    stk::mesh::Entity new_entity = new_nodes[i];

    std::vector<stk::mesh::EntitySideComponent> sides;

    //loop over the relations and check to see if they are in the closure
    for (stk::mesh::EntityRank irank = stk::topology::END_RANK;
        irank != stk::topology::BEGIN_RANK; )
    {
      --irank;

      int num_rels = mesh.num_connectivity(entity, irank);
      stk::mesh::Entity const * relations = mesh.begin(entity, irank);
      stk::mesh::ConnectivityOrdinal const *relation_ordinals = mesh.begin_ordinals(entity, irank);
      for (int j = num_rels - 1; j >= 0; --j)
      {
        stk::mesh::Entity current_entity = relations[j];
        unsigned side_ordinal = relation_ordinals[j];

        if (mesh.in_receive_ghost(mesh.entity_key(current_entity))) {
          // TODO deleteing the ghost triggers a logic error at the
          // end of the NEXT modification cycle.  We need to fix this!
          //mesh.destroy_entity(current_entity);
          continue;
        }
        // check if has relation in closure
        else if ( std::binary_search(closure.begin(),
                                     closure.end(),
                                     current_entity,
                                     lesser) )
        {
          sides.push_back(stk::mesh::EntitySideComponent(current_entity,side_ordinal));
        }
      }
    }

    //loop over the sides and break the relations between the old nodes
    //and set up the relations with the new
    for ( std::vector<stk::mesh::EntitySideComponent>::iterator
          itr = sides.begin(); itr != sides.end(); ++itr) {
      mesh.destroy_relation(itr->entity, entity, itr->side_ordinal);
      mesh.declare_relation(itr->entity, new_entity, itr->side_ordinal);
    }

    //copy non-induced part membership from nodes[i] to new_nodes[i]
      //there are NO non-induced parts for this example

    //copy field data from nodes[i] to new_nodes[i]
    mesh.copy_entity_fields( entity, new_entity);

    if (mesh.has_no_relations(entity)) {
      mesh.destroy_entity(entity);
    }

    if (mesh.has_no_relations(new_entity)) {
      mesh.destroy_entity(new_entity);
    }
  }
}

void communicate_and_create_shared_nodes( stk::mesh::BulkData & mesh,
                                          stk::mesh::EntityVector   & nodes,
                                          stk::mesh::EntityVector   & new_nodes)
{

  stk::CommAll comm(mesh.parallel());

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk::mesh::Entity node = nodes[i];
    stk::mesh::Entity new_node = new_nodes[i];

    stk::mesh::PairIterEntityComm entity_comm = mesh.entity_comm_sharing(mesh.entity_key(node));

    for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

      int proc = entity_comm.first->proc;
      comm.send_buffer(proc).pack<stk::mesh::EntityKey>(mesh.entity_key(node))
        .pack<stk::mesh::EntityKey>(mesh.entity_key(new_node));

    }
  }

  comm.allocate_buffers( mesh.parallel_size()/4 );

  for (size_t i = 0; i < nodes.size(); ++i) {
    stk::mesh::Entity node = nodes[i];
    stk::mesh::Entity new_node = new_nodes[i];

    stk::mesh::PairIterEntityComm entity_comm = mesh.entity_comm_sharing(mesh.entity_key(node));

    for (; entity_comm.first != entity_comm.second; ++entity_comm.first) {

      int proc = entity_comm.first->proc;
      comm.send_buffer(proc).pack<stk::mesh::EntityKey>(mesh.entity_key(node))
                            .pack<stk::mesh::EntityKey>(mesh.entity_key(new_node));

    }
  }

  comm.communicate();

  const stk::mesh::PartVector no_parts;

  for (int process = 0; process < mesh.parallel_size(); ++process) {
    stk::mesh::EntityKey old_key;
    stk::mesh::EntityKey new_key;

    while ( comm.recv_buffer(process).remaining()) {

      comm.recv_buffer(process).unpack<stk::mesh::EntityKey>(old_key)
                               .unpack<stk::mesh::EntityKey>(new_key);

      stk::mesh::Entity old_entity = mesh.get_entity(old_key);
      stk::mesh::Entity new_entity = mesh.declare_entity(new_key.rank(), new_key.id(), no_parts);

      nodes.push_back(old_entity);
      new_nodes.push_back(new_entity);

    }
  }
}

} // empty namespace

void separate_and_skin_mesh(
    stk::mesh::MetaData & fem_meta,
    stk::mesh::BulkData & mesh,
    stk::mesh::Part     & skin_part,
    std::vector< stk::mesh::EntityId > elements_to_separate,
    const stk::mesh::EntityRank rank_of_element
    )
{
  stk::mesh::EntityVector entities_to_separate;

  //select the entity only if the current process in the owner
  for (std::vector< stk::mesh::EntityId>::const_iterator itr = elements_to_separate.begin();
      itr != elements_to_separate.end(); ++itr)
  {
    stk::mesh::Entity element = mesh.get_entity(rank_of_element, *itr);
    if (mesh.is_valid(element) && mesh.parallel_owner_rank(element) == mesh.parallel_rank()) {
      entities_to_separate.push_back(element);
    }
  }

  stk::mesh::EntityVector entities_closure;
  stk::mesh::find_closure(mesh,
      entities_to_separate,
      entities_closure);

  stk::mesh::Selector select_owned = fem_meta.locally_owned_part();

  stk::mesh::EntityVector nodes;
  find_owned_nodes_with_relations_outside_closure(mesh, entities_closure, select_owned, nodes);

  //ask for new nodes to represent the copies
  std::vector<size_t> requests(fem_meta.entity_rank_count(), 0);
  requests[NODE_RANK] = nodes.size();

  mesh.modification_begin();

  // generate_new_entities creates new blank entities of the requested ranks
  stk::mesh::EntityVector new_nodes;
  mesh.generate_new_entities(requests, new_nodes);

  //communicate and create new nodes everywhere the old node is shared
  communicate_and_create_shared_nodes(mesh, nodes, new_nodes);

  copy_nodes_and_break_relations(mesh, entities_closure, nodes, new_nodes);

  mesh.modification_end();

  stk::mesh::PartVector add_parts(1,&skin_part);
  stk::mesh::skin_mesh(mesh, add_parts);

  return;
}
