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
#include <iostream>
#include <sstream>
#include <algorithm>

#include <mesh/Entity.hpp>
#include <mesh/EntityType.hpp>
#include <mesh/Mesh.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

bool Relation::operator < ( const Relation & r ) const
{
  entity_key_type lhs =   m_key ;
  entity_key_type rhs = r.m_key ;

  if ( lhs == rhs ) {
    lhs = NULL !=   m_entity ?   m_entity->key() : 0 ;
    rhs = NULL != r.m_entity ? r.m_entity->key() : 0 ;
  }

  return lhs < rhs ;
}

Relation::Relation( Entity & e , entity_id_type id )
  : m_key( entity_key( e.entity_type(), id) ) , m_entity(&e)
{}

//----------------------------------------------------------------------

struct LessRelation {
  inline
  bool operator()( const Relation & lhs , const Relation & rhs ) const
    { return lhs.operator<(rhs); }
};

struct LessRelationAttr {
  inline
  bool operator()( const Relation & lhs , const Relation & rhs ) const
    { return lhs.key() < rhs.key(); }

  inline
  bool operator()( const Relation & lhs , const entity_key_type rhs ) const
    { return lhs.key() < rhs ; }
};

//----------------------------------------------------------------------

std::ostream &
print_entity_key( std::ostream & os , unsigned type , entity_id_type id )
{
  const char * const name = entity_type_name( type );
  return os << name << "[" << id << "]" ;
}

std::ostream &
print_entity_key( std::ostream & os , entity_key_type key )
{
  const unsigned type     = entity_rank(key);
  const entity_id_type id = entity_id(key);
  return print_entity_key( os , type , id );
}

std::ostream &
print_relation( std::ostream & os ,
               entity_key_type patch_key ,
               entity_key_type entity_key )
{
  const unsigned      patch_type   = entity_rank( patch_key );
  const unsigned      entity_type  = entity_rank( entity_key );
  const entity_id_type local  = entity_id( patch_key );
  const entity_id_type global = entity_id( entity_key );

  const char * patch_name = entity_type_name( patch_type );

  os << patch_name << "[" << local << "->" ;

  if ( patch_type == entity_type ) {
    os << global ;
  }
  else {
    const char * type_name = entity_type_name( entity_type );
    os << type_name << "[" << global << "] " ;
  }

  os << global << "]" ;

  return os ;
}

std::ostream &
operator << ( std::ostream & os , const Relation & con )
{
  const char * name = entity_type_name( con.entity_type() );
  const unsigned local_id     = con.identifier();
  const Entity * const entity = con.entity();

  os << name << "[" << local_id << "->" ;
  if ( entity != NULL ) { os << entity->identifier(); }
  else                  { os << "?" ; }
  os << "]" ;

  return os ;
}

std::ostream &
print_entity( std::ostream & os , const std::string & lead , const Entity & e )
{
  print_entity_key( os , e.key() );
  os << " Owner(P" << e.owner_rank() << ") Relationions {" ;

  for ( RelationSpan con = e.relations() ; con ; ++con ) {
    os << std::endl << lead << "  " << *con ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

namespace {
const entity_key_type & zero_key()
{ static entity_key_type z = 0 ; return z ; }
}

Entity::Entity()
  : SetvMember< entity_key_type >( zero_key() ),
    m_relation(), m_kernel(), m_kernel_ord(0), m_owner_rank(0)
{}

Entity::~Entity()
{
  if ( m_kernel ) {
    std::ostringstream msg ;
    msg << "phdmesh::Entity::~Entity() ERROR: " ;
    print_entity_key( msg , key() );
    msg << " is still a member of a phdmesh::Kernel" ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

namespace {

RelationSpan con_span( const RelationSet & con ,
                      const entity_key_type lo_attr ,
                      const entity_key_type hi_attr )
{
  RelationSet::const_iterator i = con.begin();
  RelationSet::const_iterator e = con.end();

  i = std::lower_bound( i , e , lo_attr , LessRelationAttr() );
  e = std::lower_bound( i , e , hi_attr , LessRelationAttr() );

  return RelationSpan( i , e );
}

}

RelationSpan Entity::relations( unsigned et ) const
{
  const entity_key_type lo_key = entity_key( et , 0 );
  const entity_key_type hi_key = entity_key( et + 1 , 0 );

  return con_span( m_relation , lo_key , hi_key );
}

RelationSpan Entity::relations( unsigned et, entity_id_type id ) const
{
  const entity_key_type lo_key = entity_key( et , id );
  const entity_key_type hi_key = entity_key( et , id );

  return con_span( m_relation , lo_key , hi_key );
}

//----------------------------------------------------------------------

namespace {

void throw_required_unique( Entity & e_hi ,
                            Entity & e_lo ,
                            Entity & e_lo_exist ,
                            const entity_id_type identifier ,
                            const char * required_unique_by )
{
  static const char method_name[] = "phdmesh::Mesh::declare_relation" ;

  std::ostringstream msg ;

  msg << method_name << "( " ;
  print_entity_key( msg , e_hi.key() );
  msg << "->( " ;
  msg << identifier ;
  msg << " , " ;
  print_entity_key( msg , e_lo.key() );
  msg << " ) , " ;
  if ( required_unique_by ) {
    msg << required_unique_by ;
  }
  else {
    msg << "NULL" ;
  }
  msg << " ) FAILED : ALREADY HAS " ;
  print_entity_key( msg , e_hi.key() );
  msg << "->( " ;
  msg << identifier ;
  msg << " , " ;
  print_entity_key( msg , e_lo_exist.key() );
  msg << " )" ;

  throw std::invalid_argument(msg.str());
}

void throw_require_different( Entity & e1 , Entity & e2 ,
                              const char * required_unique_by )
{
  static const char method_name[] = "phdmesh::Mesh::declare_relation" ;

  std::ostringstream msg ;

  msg << method_name << "( " ;
  print_entity_key( msg , e1.key() );
  msg << " , " ;
  print_entity_key( msg , e2.key() );
  msg << " , " ;
  if ( required_unique_by ) {
    msg << required_unique_by ;
  }
  else {
    msg << "NULL" ;
  }
  msg << " ) FAILED : Cannot relation entities of the same type; " ;
  msg << " must introduce a relationing entity of a higher rank." ;

  throw std::invalid_argument(msg.str());
}

}

void Mesh::declare_relation( Entity & e1 , Entity & e2 ,
                               const entity_id_type identifier ,
                               const char * required_unique_by )
{
  const unsigned e1_type = e1.entity_type();
  const unsigned e2_type = e2.entity_type();

  if ( e1_type == e2_type ) {
    throw_require_different( e1 , e2 , required_unique_by );
  }

  Entity & e_hi = e1_type == e2_type ? e1 : ( e1_type < e2_type ? e2 : e1 );
  Entity & e_lo = e1_type == e2_type ? e2 : ( e1_type < e2_type ? e1 : e2 );

  {
    const RelationSet::iterator e = e_hi.m_relation.end();
          RelationSet::iterator i = e_hi.m_relation.begin();

    const Relation hi_to_lo( e_lo , identifier );

    i = std::lower_bound( i , e , hi_to_lo , LessRelation() );

    if ( e == i || hi_to_lo != *i ) { // Not a duplicate

      if ( e != i ) {
        const entity_key_type attr_hi_to_lo = hi_to_lo.key();
        const entity_key_type attr_existing = i->key();

        if ( attr_hi_to_lo == attr_existing ) {
          throw_required_unique( e_hi, e_lo, *i->entity(), identifier,
                                 required_unique_by);
        }
      }

      e_hi.m_relation.insert( i , hi_to_lo );
    }
  }

  {
    const RelationSet::iterator e = e_lo.m_relation.end();
          RelationSet::iterator i = e_lo.m_relation.begin();

    const Relation lo_to_hi( e_hi , identifier );

    i = std::lower_bound( i , e , lo_to_hi , LessRelation() );

    if ( e == i || lo_to_hi != *i ) { // Not a duplicate
      e_lo.m_relation.insert( i , lo_to_hi );
    }
  }
}

void Mesh::declare_relation( Entity & e ,
                               const std::vector<Relation> & con ,
                               const char * required_unique_by )
{
  for ( std::vector<Relation>::const_iterator
        i = con.begin() ; i != con.end() ; ++i ) {

    const unsigned r_ident  =   i->identifier();
          Entity & r_entity = * i->entity();

    declare_relation( e , r_entity , r_ident , required_unique_by );
  }
}

void Mesh::destroy_relation( Entity & e1 , Entity & e2 )
{
  RelationSet::iterator i ;

  for ( i = e1.m_relation.end() ; i != e1.m_relation.begin() ; ) {
    --i ;
    if ( & e2 == i->entity() ) { i = e1.m_relation.erase( i ); }
  }

  for ( i = e2.m_relation.end() ; i != e2.m_relation.begin() ; ) {
    --i ;
    if ( & e1 == i->entity() ) { i = e2.m_relation.erase( i ); }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

