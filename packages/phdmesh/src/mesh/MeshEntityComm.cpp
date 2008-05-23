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
#include <mesh/FieldData.hpp>
#include <mesh/EntityComm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

EntityComm::~EntityComm() {}

const char * EntityComm::name() const
{
  static const char my_name[] = "phdmesh::EntityComm" ;
  return my_name ;
}

//----------------------------------------------------------------------

void EntityComm::send_entity(
  CommBuffer & buffer ,
  const Mesh & receive_mesh ,
  const EntityProcSet::const_iterator ibeg ,
  const EntityProcSet::const_iterator iend ) const
{
  pack_entity( buffer , receive_mesh , ibeg , iend );
}

void EntityComm::receive_entity(
  CommBuffer              & buffer ,
  Mesh                    & receive_mesh ,
  const unsigned            send_source ,
  EntityProcSet & receive_info ) const
{
  entity_key_type       key ;
  unsigned              owner_rank ;
  std::vector<Part*>    parts ;
  std::vector<Relation>  relations ;
  std::vector<unsigned> send_destinations ;

  unpack_entity( buffer , receive_mesh ,
                 key , owner_rank ,
                 parts , relations , send_destinations );

  EntityProc ep ;

  {
    Entity & entity = receive_mesh.declare_entity( key );
    receive_mesh.change_entity_parts( entity , parts );
    receive_mesh.change_entity_owner( entity , owner_rank );
    ep.first = & entity ;
  }

  receive_mesh.declare_relation( *ep.first , relations , name() );

  ep.second = send_source ;

  receive_info.push_back( ep );
}

//----------------------------------------------------------------------

void EntityComm::pack_entity(
  CommBuffer & buf ,
  const Mesh & recv_mesh ,
  const EntityProcSet::const_iterator ibeg ,
  const EntityProcSet::const_iterator iend ) const
{
  const Schema      & recv_schema = recv_mesh.schema();
  const Entity      & entity      = * ibeg->first ;
  const entity_key_type key         = entity.key();
  const unsigned      owner_rank  = entity.owner_rank();
  const Kernel      & kernel      = entity.kernel();
  const Schema      & send_schema = kernel.mesh().schema();
  const bool          same_schema = & send_schema == & recv_schema ;
  const unsigned      dest_size   = std::distance( ibeg , iend );
        RelationSpan   rel         = entity.relations();

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

  buf.pack<entity_key_type>( key );
  buf.pack<unsigned>( owner_rank );

  // Parts:
  {
    const unsigned n = part_ordinals.size();
    buf.pack<unsigned>( n );
    buf.pack<unsigned>( & part_ordinals[0] , n );
  }

  // Relationships:
  {
    const unsigned rel_size = rel.size();
    buf.pack<unsigned>( rel_size );
    entity_key_type rel_data[2] ;
    for ( ; rel ; ++rel ) {
      rel_data[0] = rel->entity()->key();
      rel_data[1] = rel->key();
      buf.pack<entity_key_type>( rel_data , 2 );
    }
  }

  // Processors:
  buf.pack<unsigned>( dest_size );
  for ( EntityProcSet::const_iterator i = ibeg ; i != iend ; ++i ) {
    buf.pack<unsigned>( i->second );
  }
}

void EntityComm::unpack_entity(
  CommBuffer            & recv_buf ,
  const Mesh            & recv_mesh ,
  entity_key_type       & key ,
  unsigned              & owner_rank ,
  std::vector<Part*>    & parts ,
  std::vector<Relation>  & relations ,
  std::vector<unsigned> & send_destinations ) const
{
  const PartSet & recv_schema_parts = recv_mesh.schema().get_parts();

  parts.clear();
  relations.clear();
  send_destinations.clear();

  recv_buf.unpack<entity_key_type>( key );

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

  // Relationships:
  {
    unsigned rel_size ; recv_buf.unpack<unsigned>( rel_size );

    relations.reserve( rel_size );

    entity_key_type rel_data[2] ;

    for ( unsigned j = 0 ; j < rel_size ; ++j ) {
      recv_buf.unpack<entity_key_type>( rel_data , 2 );

      Entity * rel_entity = recv_mesh.get_entity( rel_data[0] );

      if ( rel_entity != NULL ) {
        Relation rel( *rel_entity , (unsigned) rel_data[1] );
        relations.push_back( rel );
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

void EntityComm::pack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  const Kernel & kernel = entity.kernel();
  const Mesh   & mesh   = kernel.mesh();
  const Schema & schema = mesh.schema();

  const std::vector< FieldBase * > & fields = schema.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    buf.pack<unsigned>( size );
  }

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    if ( size ) {
      unsigned char * const ptr = (unsigned char *) field_data( f , entity );
      buf.pack<unsigned char>( ptr , size );
    }
  }
}

void EntityComm::unpack_field_values(
  CommBuffer & buf , Entity & entity ) const
{
  static const char method[] = "phdmesh::EntityComm::unpack_field_values" ;

  const Kernel & kernel = entity.kernel();
  const Mesh   & mesh   = kernel.mesh();
  const Schema & schema = mesh.schema();

  const std::vector< FieldBase * > & fields = schema.get_fields();

  const std::vector< FieldBase * >::const_iterator i_end = fields.end();
  const std::vector< FieldBase * >::const_iterator i_beg = fields.begin();

  std::vector< FieldBase * >::const_iterator i ;

  std::ostringstream msg ;
  bool ok = true ;

  for ( i = i_beg ; i_end != i ; ) {
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    unsigned recv_data_size ; buf.unpack<unsigned>( recv_data_size );
    if ( size != recv_data_size ) {
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
    const FieldBase & f = **i ; ++i ;
    const unsigned size = field_data_size( f , kernel );
    if ( size ) {
      unsigned char * ptr = (unsigned char *) field_data( f , entity );
      buf.unpack<unsigned char>( ptr , size );
    }
  }
}

//----------------------------------------------------------------------

}

