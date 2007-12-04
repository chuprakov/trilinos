/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 */

#include <sstream>
#include <algorithm>
#include <stdexcept>

#include <util/ParallelReduce.hpp>
#include <util/OctTreeOps.hpp>

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>
#include <mesh/EntityComm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

namespace {

//----------------------------------------------------------------------

void global_coordinate_bounds( Mesh & M ,
                               const Field<double,1> & node_coord ,
                               double * const bounds )
{
  const Schema & schema = M.schema();
  Part & owns_part = schema.owns_part();

  double min[3] , max[3] ;
  {
    const double f_max = std::numeric_limits<float>::max();
    min[0] = min[1] = min[2] = f_max ;
    max[0] = max[1] = max[2] = - f_max ;
  }

  // Mininum and maximum for all nodes:

  const KernelSet & node_kernels = M.kernels( Node );
  const KernelSet::const_iterator i_end = node_kernels.end();
        KernelSet::const_iterator i     = node_kernels.begin();

  while ( i != i_end ) {
    const Kernel & kernel = *i ; ++i ;
    if ( kernel.has_superset( owns_part ) ) {
      const unsigned number = kernel.size();
      const double * coord = kernel.data( node_coord );

      for ( unsigned j = 0 ; j < number ; ++j ) {
        Min<3>( min , coord );
        Max<3>( max , coord );
        coord += 3 ;
      }
    }
  }

  all_reduce( M.parallel() , ReduceMin<3>( min ) &
                                  ReduceMax<3>( max ) );

  // A bounding cube:

  double center[3] , disp[3] ;
  center[0] = ( max[0] + min[0] ) * 0.5 ;
  center[1] = ( max[1] + min[1] ) * 0.5 ;
  center[2] = ( max[2] + min[2] ) * 0.5 ;
  disp[0]   = ( max[0] - min[0] ) * 0.5 ;
  disp[1]   = ( max[1] - min[1] ) * 0.5 ;
  disp[2]   = ( max[2] - min[2] ) * 0.5 ;

  double disp_max = disp[0] ;
  if ( disp_max < disp[1] ) { disp_max = disp[1] ; }
  if ( disp_max < disp[2] ) { disp_max = disp[2] ; }
  {
    const double f_eps = std::numeric_limits<float>::epsilon();
    disp_max *= 1.0 + f_eps ;
  }

  bounds[0] = center[0] - disp_max ;
  bounds[1] = center[1] - disp_max ;
  bounds[2] = center[2] - disp_max ;
  bounds[3] = disp_max * 2 ;
}

//----------------------------------------------------------------------

OctTreeKey elem_key( const double * const bounds ,
                     const Field<double,1> & node_coord ,
                     Entity & elem )
{
  enum { shift = OctTreeKey::BitsPerWord - OctTreeKey::MaxDepth };
  enum { ncell = 1 << OctTreeKey::MaxDepth };

  const double s = ncell / bounds[3] ;

  ConnectSpan elem_nodes = elem.connections( Node );

  const unsigned num_nodes = elem_nodes.size();

  double centroid[3] = { 0 , 0 , 0 };

  while ( elem_nodes ) {
    Entity & node = * elem_nodes->entity(); ++elem_nodes ;

    const double * const coord = node.data( node_coord );

    centroid[0] += coord[0] ;
    centroid[1] += coord[1] ;
    centroid[2] += coord[2] ;
  }

  unsigned oct_coord[3] ;

  oct_coord[0] = (unsigned) ( ( centroid[0] / num_nodes - bounds[0] ) * s );
  oct_coord[1] = (unsigned) ( ( centroid[1] / num_nodes - bounds[1] ) * s );
  oct_coord[2] = (unsigned) ( ( centroid[2] / num_nodes - bounds[2] ) * s );

  oct_coord[0] <<= shift ;
  oct_coord[1] <<= shift ;
  oct_coord[2] <<= shift ;

  return hsfc3d( OctTreeKey::MaxDepth , oct_coord );
}

//----------------------------------------------------------------------

void global_element_cuts( Mesh & M ,
                          const double * const bounds ,
                          const Field<double,1> & node_coord ,
                          const Field<float,1>  * const elem_weight_field ,
                                OctTreeKey      * const cut_begin )
{
  const Schema & schema = M.schema();
  const Part & owns_part = schema.owns_part();

  const KernelSet & elem_kernels = M.kernels( Element );
  const KernelSet::const_iterator i_end = elem_kernels.end();
        KernelSet::const_iterator i     = elem_kernels.begin();

  unsigned number = 0 ;

  while ( i != i_end ) {
    const Kernel & kernel = *i ; ++i ;
    if ( kernel.has_superset( owns_part ) ) {
      number += kernel.size();
    }
  }

  std::vector<OctTreeKey> elem_keys( number );
  std::vector<float>      elem_weights( elem_weight_field ? number : 0 );

  i = elem_kernels.begin();

  unsigned count = 0 ;

  while ( i != i_end ) {
    const Kernel & kernel = *i ; ++i ;
    if ( kernel.has_superset( owns_part ) ) {

      const Kernel::iterator j_end = kernel.end();
            Kernel::iterator j     = kernel.begin();

      while ( j != j_end ) {
        Entity & elem = **j ; ++j ;

        elem_keys[count] = elem_key( bounds , node_coord , elem );

        if ( elem_weight_field ) {
          const float one = 1.0 ;
          const float * const w = elem.data( *elem_weight_field );
          elem_weights[count] = w ? *w : one ;
        }

        ++count ;
      }
    }
  }

  // Generate partitioning cuts for the global set of keys

  const OctTreeKey * const keys = & elem_keys[0] ; 

  const float * const weights =
    elem_weight_field ? & elem_weights[0] : (float *) NULL ;

  oct_tree_partition_fine( M.parallel() , number , keys , weights ,
                           cut_begin );
}

//----------------------------------------------------------------------
// Check each non-aura "other" entity against the attached uses-entities.
// Each entity will be rebalanced to every processor that it 'uses'.

void rebal_other_entities( Mesh & M , EntityProcSet & rebal )
{
  static const char method[] = "phdmesh::rebal_other_entities" ;

  sort_unique( rebal );

  const Schema & schema = M.schema();
  const Part & uses_part = schema.uses_part();

  const KernelSet::iterator k_end = M.kernels(Other).end();
        KernelSet::iterator k     = M.kernels(Other).begin();

  std::vector<unsigned> rebal_procs ;
  EntityProcSet rebal_tmp ;

  while ( k != k_end ) {
    Kernel & kernel = *k ; ++k ;

    if ( kernel.has_superset( uses_part ) ) {

      const Kernel::iterator i_end = kernel.end();
            Kernel::iterator i     = kernel.begin();

      while ( i != i_end ) {
        Entity * const entity = *i ; ++i ;

        rebal_procs.clear();

        for ( ConnectSpan j = entity->connections(); j ; ++j ) {
          if ( j->entity_type() == Other ) {
            std::ostringstream msg ;
            msg << "P" << M.parallel_rank();
            msg << ": " << method ;
            msg << " FAILED, Cannot rebalance " ;
            print_entity_key( msg , entity->key() );
            msg << "->{ " ;
            msg << *j ;
            msg << " }" ;
            throw std::runtime_error( msg.str() );
          }
        }

        for ( ConnectSpan j = entity->connections(); j ; ++j ) {
          Entity * const con_entity = j->entity();
          EntityProcSet::const_iterator ir =
            lower_bound( rebal , *con_entity );
          while ( ir != rebal.end() && ir->first == con_entity ) {
            rebal_procs.push_back( ir->second );
          }
        }

        std::vector<unsigned>::iterator ip = rebal_procs.begin();
        std::vector<unsigned>::iterator ep = rebal_procs.end();

        std::sort( ip , ep );
        ip = std::unique( ip , ep );
        rebal_procs.erase( ip , ep );

        for ( ip = rebal_procs.begin() ; ip != rebal_procs.end() ; ++ip ) {
          EntityProc tmp( entity , *ip );
          rebal_tmp.push_back( tmp );
          for ( ConnectSpan j = entity->connections(); j ; ++j ) {
            Entity * const con_entity = j->entity();
            EntityProc con_tmp( con_entity , *ip );
            rebal_tmp.push_back( con_tmp );
          }
        }
      }
    }
  }

  rebal.insert( rebal.end() , rebal_tmp.begin() , rebal_tmp.end() );

  sort_unique( rebal );
}

// Check each non-aura element-entity against the attached elements:

void rebal_elem_entities( Mesh & M ,
                          EntityType t ,
                          EntityProcSet & rebal )
{
  const Schema & schema = M.schema();
  const Part & uses_part = schema.uses_part();

  const EntitySet::iterator i_end = M.entities(t).end();
        EntitySet::iterator i     = M.entities(t).begin();

  std::vector<unsigned> rebal_procs ;

  while ( i != i_end ) {
    Entity * const entity = & *i ; ++i ;

    if ( entity->kernel().has_superset( uses_part ) ) {

      rebal_procs.clear();

      ConnectSpan elements = entity->connections( Element );

      while ( elements ) {
        Entity & elem = * elements->entity(); ++elements ;
        const unsigned p = elem.owner_rank();
        rebal_procs.push_back( p );
      }

      std::vector<unsigned>::iterator j = rebal_procs.begin();
      std::vector<unsigned>::iterator e = rebal_procs.end();

      std::sort( j , e );
      j = std::unique( j , e );
      rebal_procs.erase( j , e );

      for ( j = rebal_procs.begin() ; j != rebal_procs.end() ; ++j ) {
        EntityProc tmp( entity , *j );
        rebal.push_back( tmp );
      }

      M.change_entity_owner( *entity , rebal_procs.back() );
    }
  }
}

class RebalanceComm : public EntityComm {
public:
  ~RebalanceComm() {}
  RebalanceComm() {}

  const char * name() const ;

  void send_entity( CommBuffer & , const Mesh & ,
                    const EntityProcSet::const_iterator ,
                    const EntityProcSet::const_iterator ) const ;

  void receive_entity(
    CommBuffer              & buffer ,
    Mesh                    & receive_mesh ,
    const unsigned            send_source ,
    EntityProcSet & receive_info ) const ;

private:
  RebalanceComm( const RebalanceComm & );
  RebalanceComm & operator = ( const RebalanceComm & );
};

const char * RebalanceComm::name() const
{ static const char n[] = "phdmesh::RebalanceComm" ; return n ; }

void RebalanceComm::send_entity(
  CommBuffer & buf ,
  const Mesh & recv_mesh ,
  const EntityProcSet::const_iterator ib ,
  const EntityProcSet::const_iterator ie ) const 
{
  // Only the current owner sends an entity:

  Entity & entity = * ib->first ;
  Part & owns_part = recv_mesh.schema().owns_part();

  if ( entity.kernel().has_superset( owns_part ) ) {
    pack_entity(       buf , recv_mesh , ib , ie );
    pack_field_values( buf , entity );
  }
}

void RebalanceComm::receive_entity(
  CommBuffer              & buffer ,
  Mesh                    & receive_mesh ,
  const unsigned            send_source ,
  EntityProcSet & receive_info ) const
{
  EntityType            entity_type ;
  unsigned long         entity_id ;
  unsigned              owner_rank ;
  std::vector<Part*>    parts ;
  std::vector<Connect>  connections ;
  std::vector<unsigned> send_dest ;

  unpack_entity( buffer , receive_mesh ,
                 entity_type , entity_id , owner_rank ,
                 parts , connections , send_dest );

  EntityProc ep ;

  ep.first = & receive_mesh.declare_entity( entity_type , entity_id ,
                                            parts , owner_rank );

  receive_mesh.declare_connection( *ep.first , connections , name() );

  ep.second = send_source ;

  std::vector<unsigned>::iterator ip ;
  for ( ip = send_dest.begin() ; ip != send_dest.end() ; ++ip ) {
    ep.second = *ip ;
    receive_info.push_back( ep );
  }

  unpack_field_values( buffer , * ep.first );
}

void destroy_not_retained( Mesh & M , EntityProcSet & rebal )
{
  const unsigned p_rank = M.parallel_rank();

  // Iterating backwards, the 'rebal' array must be sorted and unique

  EntityProcSet::iterator i = rebal.end();

  while ( i != rebal.begin() ) {
    const EntityProcSet::iterator j = i ;

    Entity * const entity = (--i)->first ;

    bool remove = true ;
    bool flag ;

    do {
      if ( p_rank == i->second ) { remove = false ; }

      flag = i != rebal.begin();

      if ( flag ) {
        EntityProcSet::iterator k = i ;

        flag = entity == (--k)->first ;

        if ( flag ) { --i ; }
      }

    } while ( flag );

    if ( remove ) {
      i = rebal.erase( i , j );

      M.destroy_entity( entity );
    }
  }
}

//----------------------------------------------------------------------

void remove_aura( Mesh & M )
{
  Part & uses_part = M.schema().uses_part();

  for ( Span< std::vector< EntityProc >::const_iterator >
        span( M.aura_range() ) ; span ; ++span ) {
    Entity * const e = span->first ;
    if ( e->kernel().has_superset( uses_part ) ) {
      std::string msg("phdmesh::remove_aura corrupted aura");
      throw std::runtime_error(msg);
    }
    M.destroy_entity( e );
  }

  const EntityProcSet tmp ;
  M.set_aura( tmp , tmp );
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Dynamic load balancing algorithm:
//
//  [1] Determine the new load balance partitioning for the
//      global set of locally owned elements (avoid duplicatation).
//  [2] Set all element's owner_rank to the rebalance processor.
//  [3] Each non-aura node, edge, and face:
//      - Should be completely surrounded by its elements due
//        to aura elements.
//      - If attached to an owned element then add it to the rebalance list.
//      - If not attached to an element that will remain on this
//        processor then add it to the deletion list.
//      - If attached to an element that will remain on this
//        processor then add to the resident list with all
//        using elements' processors.
//  [4] Delete all aura entities, are done with them
//  [5] Add owned elements with a remote owner_rank to the relocation list.
//      Add elements with a remote owner_rank to the deletion list.
//  [6] Pack relocation list (forward order).
//      - Include processor ranks for all elements using the entities.
//      - Include all field values.
//  [7] Delete entities from the deletion list (reverse order).
//  [8] Communicate relocation and unpack.
//      - Add to the resident list with all incoming using elements'
//        processors.
//  [9] Determine ownership and sharing by processing the resident list.
//  [10] Create new aura layer.
//
//----------------------------------------------------------------------

void comm_mesh_rebalance( Mesh & M ,
                          const Field<double,1> & node_coord_field ,
                          const Field<float,1>  * const elem_weight_field ,
                          std::vector<OctTreeKey> & cut_keys )
{
  const Schema & schema  = M.schema();
  Part * const uses_part = & schema.uses_part();
  Part * const owns_part = & schema.owns_part();

  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  //--------------------------------------------------------------------
  // The node_coord_field must be up to date on all processors
  // so that the element oct tree keys are parallel consistent.
  // It is assumed that the shared node_coord_field values are
  // already consistent.
  {
    const Field<void,0> * const ptr = & node_coord_field ;
    std::vector< const Field<void,0> *> tmp ;
    tmp.push_back( ptr );
    const EntityProcSet & aura_domain = M.aura_domain();
    const EntityProcSet & aura_range  = M.aura_range();
    comm_mesh_field_values( M , aura_domain , aura_range , tmp , false );
  }
  //--------------------------------------------------------------------
  // Generate global oct-tree keys for local element centroids
  // and cuts for the global element centroids.

  double bounds[4] ;

  global_coordinate_bounds( M , node_coord_field , bounds );

  cut_keys.assign( p_size , OctTreeKey() );

  OctTreeKey * const cut_begin = & cut_keys[0] ;
  OctTreeKey * const cut_first = cut_begin + 1 ;
  OctTreeKey * const cut_end   = cut_begin + p_size ;

  global_element_cuts( M , bounds , node_coord_field ,
                                    elem_weight_field , cut_begin );

  //--------------------------------------------------------------------
  // Mapping of *all* elements to load balanced processor,
  // even the aura elements.
  // This requires that the node coordinates on the aura
  // elements be up to date.

  {
    std::vector< const Field<void,0> * > tmp ;
    const Field<void,0> * const tmp_coord = & node_coord_field ;
    tmp.push_back( tmp_coord );

    const EntityProcSet & d = M.aura_domain();
    const EntityProcSet & r = M.aura_range();

    comm_mesh_field_values( M , d , r , tmp , false );
  }

  {
    const EntitySet & elem_set = M.entities( Element );
    const EntitySet::iterator i_end = elem_set.end();
          EntitySet::iterator i     = elem_set.begin();
    while ( i != i_end ) {
      Entity & elem = *i ; ++i ;

      const OctTreeKey k = elem_key( bounds , node_coord_field , elem );

      const unsigned p = std::upper_bound(cut_first, cut_end, k) - cut_first ;

      M.change_entity_owner( elem , p );
    }
  }
  //--------------------------------------------------------------------
  // Fill 'rebal' with all uses entities' rebalancing processors

  EntityProcSet rebal ;

  rebal_elem_entities( M , Node , rebal );
  rebal_elem_entities( M , Edge , rebal );
  rebal_elem_entities( M , Face , rebal );

  {
    const Part & part_uses = * uses_part ;

    const KernelSet & elem_kernels = M.kernels( Element );
    const KernelSet::const_iterator i_end = elem_kernels.end();
          KernelSet::const_iterator i     = elem_kernels.begin();

    while ( i != i_end ) {
      const Kernel & kernel = *i ; ++i ;

      if ( kernel.has_superset( part_uses ) ) {

        const Kernel::iterator j_end = kernel.end();
              Kernel::iterator j     = kernel.begin();

        while ( j != j_end ) {
          Entity * const entity = *j ; ++j ;
          const unsigned p = entity->owner_rank();
          EntityProc tmp( entity , p );
          rebal.push_back( tmp );
        }
      }
    }
  }

  // The 'other' entities rebalance based upon the entities
  // that they use.  This may lead to more sharing entities.
  // Thus 'rebal' is input and then updated.

  rebal_other_entities( M , rebal );

  // 'rebal' now contains the rebalancing (entity,processor) pairs
  // for every non-aura entity.  Can now delete the aura entities.

  remove_aura( M );

  // Copy entities to new processors according to 'rebal'.
  // Only send the owned entities.
  // Include all processors associated with the entity in 'rebal'.
  // Unpack all nodes, then all edges, then all faces, then all elements, 
  // from each processor.
  // The owner of a shared entity is the max-rank processor.
  // Add received entities to shared if more than one processor.

  {
    const RebalanceComm manager ;
    EntityProcSet recv_rebal ;

    comm_mesh_entities( manager , M , M , rebal , recv_rebal , false );

    // Destroy not-retained entities, they have been packed.
    // Remove the corresponding entries in 'rebal'

    destroy_not_retained( M , rebal );

    rebal.insert( rebal.end() , recv_rebal.begin() , recv_rebal.end() );

    sort_unique( rebal );
  }

  // The 'rebal' should contain a reference to every non-aura entity
  // on the local processor.  These references include every
  // processor on which the entity now resides, including the
  // local processor.

  { // Set parallel ownership and sharing parts.

    EntityProcSet::iterator ish ;

    for ( ish = rebal.begin() ; ish != rebal.end() ; ) {
      Entity & e = * ish->first ;

      for ( ; ish != rebal.end() && ish->first == & e ; ++ish );

      const bool is_owned = p_rank == e.owner_rank() ;

      // Change ownership.

      std::vector<Part*> add_parts ;
      std::vector<Part*> remove_parts ;

      if ( is_owned ) { add_parts.push_back( owns_part ); }
      else            { remove_parts.push_back( owns_part ); }

      M.change_entity_parts( e , add_parts , remove_parts );
    }

    // Remove references to the local processor,
    // the remaining entries define the sharing.

    for ( ish = rebal.end() ; ish != rebal.begin() ; ) {
      --ish ;
      if ( p_rank == ish->second ) { ish = rebal.erase( ish ); }
    }

    M.set_shares( rebal );
  }

  // Establish new aura

  comm_mesh_regenerate_aura( M );
}

}

