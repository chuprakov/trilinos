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

#include <stdexcept>
#include <sstream>

#include <mesh/Schema.hpp>
#include <mesh/Part.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Entity.hpp>
#include <element/LocalTopology.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const LocalTopology * get_part_local_topology( Part & p )
{
  const LocalTopology * top = NULL ;

  CSet::Span<LocalTopology> span = p.attribute<LocalTopology>();

  if ( 1 == span.size() ) {
    top = & *span ;
  }
  else if ( 1 < span.size() ) {
    std::ostringstream msg ;
    msg << "phdmesh::get_part_local_topology( " ;
    msg << p.name();
    msg << " ) ERROR, too many topologies { " ;
    for ( ; span ; ++span ) {
      msg << span->name << " " ;
    }
    msg << "}" ;
    throw std::runtime_error( msg.str() );
  }
  return top ;
}

void set_part_local_topology( Part & p , const LocalTopology * singleton )
{
  static const char method[] = "phdmesh::set_part_local_topology" ;

  CSet::Span<LocalTopology> span = p.attribute<LocalTopology>();

  const bool error_null   = singleton == NULL ;
  const bool error_size   = 1 < span.size();
  const bool error_change = span.size() == 1 && singleton != & *span ;

  if ( error_size || error_change || error_null ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( " ;
    msg << p.name();
    msg << " , " ;
    if ( error_null ) { msg << "NULL" ; }
    else              { msg << singleton->name ; }
    msg << " ) ERROR, current topologies { " ;
    for ( ; span ; ++span ) {
      msg << span->name << " " ;
    }
    msg << "}" ;
    throw std::runtime_error( msg.str() );
  }

  if ( span.empty() ) {
    p.schema().declare_part_attribute( p , singleton , false );
  }
}

//----------------------------------------------------------------------

const LocalTopology * get_local_topology( Kernel & kernel )
{
  const LocalTopology * top = NULL ;
  PartSet parts ;
  kernel.supersets( parts );

  PartSet::iterator i = parts.begin() ;

  for ( ; NULL == top && i != parts.end() ; ++i ) {
    top = get_part_local_topology( **i );
  }

  bool ok = true ;

  for ( ; ok && i != parts.end() ; ++i ) {
    const LocalTopology * const tmp = get_part_local_topology( **i );
    ok = tmp == NULL || tmp == top ;
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "phdmesh::get_local_topology( Kernel[" ;
    for ( i = parts.begin() ; i != parts.end() ; ++i ) {
      const LocalTopology * const tmp = get_part_local_topology( **i );
      msg << " " << (*i)->name();
      if ( top ) { msg << "->" << tmp->name ; }
      msg << " ] ) FAILED WITH MULTIPLE LOCAL TOPOLOGIES" ;
      throw std::runtime_error( msg.str() );
    }
  }

  return top ;
}

const LocalTopology * get_local_topology( Entity & entity )
{ return get_local_topology( entity.kernel() ); }

//----------------------------------------------------------------------

Entity & declare_element( Mesh & mesh ,
                          const LocalTopology & top ,
                          const unsigned elem_id ,
                          const unsigned node_id[] )
{
  const EntityType type =
    ! top.is_boundary ? Element : EntityType( top.topological_rank );

  const entity_key_type key = entity_key( type , elem_id );

  Entity & elem = mesh.declare_entity( key );

  for ( unsigned i = 0 ; i < top.number_node ; ++i ) {
    Entity & node = mesh.declare_entity( entity_key( Node , node_id[i] ) );
    mesh.declare_relation( elem , node , i );
  }
  return elem ;
}

Entity & declare_element( Mesh & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          const unsigned node_id[] )
{
  static const char method[] = "phdmesh::declare_element" ;

  const LocalTopology * const top = get_part_local_topology( part );

  if ( top == NULL ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "( mesh , " ;
    msg << part.name();
    msg << " , " ;
    msg << elem_id ;
    msg << " , node_id[] ) ERROR, Part does not have a local topology" ;
    throw std::runtime_error( msg.str() );
  }

  Entity & elem = declare_element( mesh , *top , elem_id , node_id );

  {
    PartSet add ;
    Part * const tmp = & part ;
    add.push_back( tmp );

    mesh.change_entity_parts( elem , add );
  }

  return elem ;
}

//----------------------------------------------------------------------

Entity & declare_element_side( Mesh   & mesh , const unsigned global_side_id ,
                               Entity & elem , const unsigned local_side_id )
{
  static const char method[] = "phdmesh::declare_element_side" ;

  const LocalTopology * const elem_top = get_local_topology( elem );

  const LocalTopology * const side_top =
    ( elem_top && local_side_id < elem_top->number_side )
    ? elem_top->side[ local_side_id ].topology : NULL ;

  if ( NULL == side_top ) {
     std::ostringstream msg ;
     msg << method << "( mesh , "
         << global_side_id
         << " , " ;
     print_entity_key( msg , elem.key() );
     msg << " , "
         << local_side_id
         << " ) FAILED" ;
     if ( NULL == elem_top ) {
       msg << " Cannot discern element topology" ;
     }
     else {
       msg << " Local side id exceeds " ;
       msg << elem_top->name ;
       msg << ".number_side = " ;
       msg << elem_top->number_side ;
     }
     throw std::runtime_error( msg.str() );
   }

  const unsigned * const side_node_map =
    elem_top->side[ local_side_id ].node ;

  EntityType side_type = (EntityType) side_top->topological_rank ;

  Entity & side =
    mesh.declare_entity( entity_key( side_type, global_side_id ) );

  RelationSpan rel = elem.relations( Node );

  for ( unsigned i = 0 ; i < side_top->number_node ; ++i ) {
    Entity & node = * rel[ side_node_map[i] ].entity();
    mesh.declare_relation( side , node , i );
  }

  mesh.declare_relation( elem , side , local_side_id );

  return side ;
}

}

