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

#include <iterator>
#include <stdexcept>
#include <sstream>

#include <mesh/Mesh.hpp>
#include <mesh/Schema.hpp>
#include <mesh/EntityManager.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

EntityManager::~EntityManager() {}

const char * EntityManager::name() const
{
  static const char my_name[] = "phdmesh::EntityManager" ;
  return my_name ;
}

Entity * EntityManager::declare_entity(
  Mesh & recv_mesh ,
  EntityType entity_type ,
  unsigned long entity_id ,
  unsigned owner_rank ,
  std::vector<Part*> & parts ,
  std::vector<Connect> & rels ) const
{
  Entity & e = recv_mesh.declare_entity( entity_type , entity_id ,
                                         parts , owner_rank );

  for ( std::vector<Connect>::iterator
        i = rels.begin() ; i != rels.end() ; ++i ) {

    const ConnectType r_type   = i->type();
    const unsigned    r_ident  = i->identifier();
          Entity    & r_entity = * i->entity();

    switch( r_type ) {
    case Uses :
      recv_mesh.declare_connection( e , r_entity , r_ident );
      break ;
    case UsedBy :
      recv_mesh.declare_connection( r_entity , e , r_ident );
      break ;
    case Anonymous :
      recv_mesh.declare_connection_anon( e , r_entity , r_ident );
      break ;
    }
  }

  return & e ;
}

//----------------------------------------------------------------------

void EntityManager::send_entity(
  CommBuffer & buffer ,
  const Mesh & receive_mesh ,
  const std::vector<EntityProc>::const_iterator ibeg ,
  const std::vector<EntityProc>::const_iterator iend ) const
{
  pack_entity( buffer , receive_mesh , ibeg , iend );
}

void EntityManager::receive_entity(
  CommBuffer              & buffer ,
  Mesh                    & receive_mesh ,
  const unsigned            send_source ,
  std::vector<EntityProc> & receive_info ) const
{
  EntityType            entity_type ;
  unsigned long         entity_id ;
  unsigned              owner_rank ;
  std::vector<Part*>    parts ;
  std::vector<Connect>  connections ;
  std::vector<unsigned> send_destinations ;

  unpack_entity( buffer , receive_mesh ,
                 entity_type , entity_id , owner_rank ,
                 parts , connections , send_destinations );

  EntityProc ep ;

  ep.first = declare_entity( receive_mesh , entity_type , entity_id ,
                             owner_rank , parts , connections );
  ep.second = send_source ;

  receive_info.push_back( ep );
}

//----------------------------------------------------------------------

void EntityManager::pack_entity(
  CommBuffer & buf ,
  const Mesh & recv_mesh ,
  const std::vector<EntityProc>::const_iterator ibeg ,
  const std::vector<EntityProc>::const_iterator iend ) const
{
  const Schema      & recv_schema = recv_mesh.schema();
  const Entity      & entity      = * ibeg->first ;
  const unsigned long key         = entity.key();
  const unsigned      owner_rank  = entity.owner_rank();
  const Kernel      & kernel      = entity.kernel();
  const Schema      & send_schema = kernel.mesh().schema();
  const bool          same_schema = & send_schema == & recv_schema ;
  const ConnectSpan   rel         = entity.connections();
  const unsigned      rel_size    = std::distance( rel.first , rel.second );
  const unsigned      dest_size   = std::distance( ibeg , iend );

  std::vector<unsigned> part_ordinals ;

  entity.kernel().supersets( part_ordinals );

  if ( ! same_schema ) { // Map the parts by name
    const PartSet & send_parts = send_schema.get_parts();
    std::vector<unsigned>::iterator i = part_ordinals.begin();

    while ( i != part_ordinals.end() ) {
      Part * const send_p = send_parts[ *i ];
      Part * const recv_p = recv_schema.get_part( send_p->name() );
      if ( recv_p != NULL ) {
        *i = recv_p->schema_ordinal(); ++i ;
      }
      else {
        i = part_ordinals.erase(i);
      }
    }
  }

  buf.pack<unsigned long>( key );
  buf.pack<unsigned>( owner_rank );

  // Parts:
  {
    const unsigned n = part_ordinals.size();
    buf.pack<unsigned>( n );
    buf.pack<unsigned>( & part_ordinals[0] , n );
  }

  // Connections:
  {
    buf.pack<unsigned>( rel_size );
    unsigned long rel_data[2] ;
    for ( std::vector<Connect>::const_iterator j =  rel.first ;
                                               j != rel.second ; ++j ) {
      rel_data[0] = j->entity()->key();
      rel_data[1] = j->attribute();
      buf.pack<unsigned long>( rel_data , 2 );
    }
  }

  // Processors:
  buf.pack<unsigned>( dest_size );
  for ( std::vector<EntityProc>::const_iterator i = ibeg ; i != iend ; ++i ) {
    buf.pack<unsigned>( i->second );
  }
}

void EntityManager::unpack_entity(
  CommBuffer            & recv_buf ,
  const Mesh            & recv_mesh ,
  EntityType            & entity_type ,
  unsigned long         & entity_id ,
  unsigned              & owner_rank ,
  std::vector<Part*>    & parts ,
  std::vector<Connect>  & connections ,
  std::vector<unsigned> & send_destinations ) const
{
  const PartSet & recv_schema_parts = recv_mesh.schema().get_parts();

  parts.clear();
  connections.clear();
  send_destinations.clear();

  {
    unsigned long key ; recv_buf.unpack<unsigned long>( key );

    entity_type = Entity::key_entity_type( key );
    entity_id   = Entity::key_identifier( key );
  }

  recv_buf.unpack<unsigned>( owner_rank );

  // Parts:
  {
    unsigned n ; recv_buf.unpack<unsigned>( n );
    for ( ; n ; --n ) {
      unsigned ordinal ; recv_buf.unpack<unsigned>( ordinal );
      Part * const p = recv_schema_parts[ ordinal ];
      parts.push_back( p );
    }
  }

  // Connections:
  {
    unsigned rel_size ; recv_buf.unpack<unsigned>( rel_size );

    connections.reserve( rel_size );

    unsigned long rel_data[2] ;

    for ( unsigned j = 0 ; j < rel_size ; ++j ) {
      recv_buf.unpack<unsigned long>( rel_data , 2 );

      Entity * rel_entity = recv_mesh.get_entity( rel_data[0] );

      if ( rel_entity != NULL ) {
        Connect rel( *rel_entity , (unsigned) rel_data[1] );
        connections.push_back( rel );
      }
    }
  }

  // Processors:
  {
    unsigned dest_size ; recv_buf.unpack<unsigned>( dest_size );

    const unsigned u_zero = 0 ;
    send_destinations.assign( dest_size , u_zero );

    unsigned * const tmp = & send_destinations[0] ;
    recv_buf.unpack<unsigned>( tmp , dest_size );
  }
}

//----------------------------------------------------------------------

void EntityManager::pack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  const EntityType entity_type = entity.entity_type();
  const Kernel & kernel = entity.kernel();
  const Mesh   & mesh   = kernel.mesh();
  const Schema & schema = mesh.schema();

  const std::vector< Field<void,0> * > & fields =
    schema.get_fields( entity_type );

  const std::vector< Field<void,0> * >::const_iterator i_end = fields.end();
  const std::vector< Field<void,0> * >::const_iterator i_beg = fields.begin();

  std::vector< Field<void,0> * >::const_iterator i ;

  for ( i = i_beg ; i_end != i ; ) {
    const Field<void,0> & f = **i ; ++i ;
    const unsigned data_size = kernel.data_size( f );
    buf.pack<unsigned>( data_size );
  }

  for ( i = i_beg ; i_end != i ; ) {
    const Field<void,0> & f = **i ; ++i ;
    const unsigned data_size = kernel.data_size( f );
    if ( data_size ) {
      unsigned char * const data = (unsigned char *) entity.data( f );
      buf.pack<unsigned char>( data , data_size );
    }
  }
}

void EntityManager::unpack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  static const char method[] = "phdmesh::EntityManager::unpack_field_values" ;

  const EntityType entity_type = entity.entity_type();
  const Kernel & kernel = entity.kernel();
  const Mesh   & mesh   = kernel.mesh();
  const Schema & schema = mesh.schema();

  const std::vector< Field<void,0> * > & fields =
    schema.get_fields( entity_type );

  const std::vector< Field<void,0> * >::const_iterator i_end = fields.end();
  const std::vector< Field<void,0> * >::const_iterator i_beg = fields.begin();

  std::vector< Field<void,0> * >::const_iterator i ;

  std::ostringstream msg ;
  bool ok = true ;

  for ( i = i_beg ; i_end != i ; ) {
    const Field<void,0> & f = **i ; ++i ;
    const unsigned data_size = kernel.data_size( f );
    unsigned recv_data_size ; buf.unpack<unsigned>( recv_data_size );
    if ( data_size != recv_data_size ) {
      if ( ok ) {
        msg << "P" << mesh.parallel_rank();
        msg << ": " << method ;
        msg << "( " ;
        print_entity_key( msg , entity.key() );
        msg << " ) FAILED, incompatible size for field {" ;
      }
      msg << " " << (*i)->name();
      ok = false ;
    }
  }

  if ( ! ok ) {
    msg << " }" ;
    throw std::runtime_error( msg.str() );
  }

  for ( i = i_beg ; i_end != i ; ) {
    const Field<void,0> & f = **i ; ++i ;
    const unsigned data_size = kernel.data_size( f );
    if ( data_size ) {
      unsigned char * data = (unsigned char *) entity.data( f );
      buf.unpack<unsigned char>( data , data_size );
    }
  }
}

//----------------------------------------------------------------------

}

