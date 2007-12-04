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
#include <mesh/Mesh.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

bool Connect::operator < ( const Connect & r ) const
{
  unsigned long lhs =   m_attr ;
  unsigned long rhs = r.m_attr ;

  if ( lhs == rhs ) {
    lhs = NULL !=   m_entity ?   m_entity->key() : 0 ;
    rhs = NULL != r.m_entity ? r.m_entity->key() : 0 ;
  }

  return lhs < rhs ;
}

Connect::Connect( Entity & e , unsigned id )
  : m_attr( attribute( e.entity_type(), id) ) , m_entity(&e)
{}

//----------------------------------------------------------------------

unsigned long Entity::create_key( EntityType type , unsigned long id )
{
  static const char method[] = "phdmesh::Entity::create_key" ;

  unsigned long key = id & KeyIdentifierMask ;

  if ( 0 == id || id != key ) {
    std::ostringstream msg ;
    msg << method << " FAILED, Identifier " << id ;
    msg << " not in range [1.." ;
    msg << ((unsigned long) KeyIdentifierMask ) << "]" ;
    throw std::invalid_argument( msg.str() );
  }

  key |= ((unsigned long)type) << KeyIdentifierDigits ;

  return key ;
}

//----------------------------------------------------------------------

struct LessConnect {
  inline
  bool operator()( const Connect & lhs , const Connect & rhs ) const
    { return lhs.operator<(rhs); }
};

struct LessConnectAttr {
  inline
  bool operator()( const Connect & lhs , const Connect & rhs ) const
    { return lhs.attribute() < rhs.attribute(); }

  inline
  bool operator()( const Connect & lhs , const unsigned rhs ) const
    { return lhs.attribute() < rhs ; }
};

//----------------------------------------------------------------------

std::ostream &
print_entity_key( std::ostream & os , EntityType type , unsigned long id )
{
  const char * const name = entity_type_name( type );
  return os << name << "[" << id << "]" ;
}

std::ostream &
print_entity_key( std::ostream & os , unsigned long key )
{
  const EntityType type  = Entity::key_entity_type( key );
  const unsigned long id = Entity::key_identifier( key );
  return print_entity_key( os , type , id );
}

std::ostream &
print_connect( std::ostream & os ,
               unsigned connect_attr ,
               unsigned long entity_key )
{
  const EntityType    type   = Connect::entity_type( connect_attr );
  const unsigned      local  = Connect::identifier( connect_attr );
  const unsigned long global = Entity::key_identifier( entity_key );

  const char * name = entity_type_name( type );

  os << name << "[" << local << "->" << global << "]" ;

  return os ;
}

std::ostream &
operator << ( std::ostream & os , const Connect & con )
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
  os << " Owner(P" << e.owner_rank() << ") Connections {" ;

  for ( ConnectSpan con = e.connections() ; con ; ++con ) {
    os << std::endl << lead << "  " << *con ;
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

namespace {
const unsigned long & zero_ul() { static unsigned long z = 0 ; return z ; }
}

Entity::Entity()
  : SetvMember<unsigned long>( zero_ul() ),
    m_connect(), m_kernel(), m_kernel_ord(0), m_owner_rank(0)
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

ConnectSpan con_span( const ConnectSet & con ,
                      const unsigned lo_attr ,
                      const unsigned hi_attr )
{
  ConnectSet::const_iterator i = con.begin();
  ConnectSet::const_iterator e = con.end();

  i = std::lower_bound( i , e , lo_attr , LessConnectAttr() );
  e = std::lower_bound( i , e , hi_attr , LessConnectAttr() );

  return ConnectSpan( i , e );
}

}

ConnectSpan Entity::connections( EntityType et ) const
{
  const unsigned lo_key = Connect::attribute( et ,   0 );
  const unsigned hi_key = Connect::attribute( et+1 , 0 );

  return con_span( m_connect , lo_key , hi_key );
}

ConnectSpan Entity::connections( EntityType et, unsigned id ) const
{
  const unsigned lo_key = Connect::attribute( et , id );
  const unsigned hi_key = Connect::attribute( et , id );

  return con_span( m_connect , lo_key , hi_key );
}

//----------------------------------------------------------------------

namespace {

void throw_required_unique( Entity & e_hi ,
                            Entity & e_lo ,
                            Entity & e_lo_exist ,
                            const unsigned identifier ,
                            const char * required_unique_by )
{
  static const char method_name[] = "phdmesh::Mesh::declare_connection" ;

  std::ostringstream msg ;

  msg << method_name << "( " ;
  print_entity_key( msg , e_hi.key() );
  msg << "->[" ;
  msg << identifier ;
  msg << "]->" ;
  print_entity_key( msg , e_lo.key() );
  msg << " , " ;
  msg << required_unique_by ;
  msg << " ) FAILED : ALREADY HAS " ;
  print_entity_key( msg , e_hi.key() );
  msg << "->[" ;
  msg << identifier ;
  msg << "]->" ;
  print_entity_key( msg , e_lo_exist.key() );

  throw std::invalid_argument(msg.str());
}

}

void Mesh::declare_connection( Entity & e1 , Entity & e2 ,
                               const unsigned identifier ,
                               const char * required_unique_by )
{
  const EntityType e1_type = e1.entity_type();
  const EntityType e2_type = e2.entity_type();

  Entity & e_hi = e1_type == e2_type ? e1 : ( e1_type < e2_type ? e2 : e1 );
  Entity & e_lo = e1_type == e2_type ? e2 : ( e1_type < e2_type ? e1 : e2 );

  {
    const ConnectSet::iterator e = e_hi.m_connect.end();
          ConnectSet::iterator i = e_hi.m_connect.begin();

    const Connect hi_to_lo( e_lo , identifier );

    i = std::lower_bound( i , e , hi_to_lo , LessConnect() );

    if ( e == i || hi_to_lo != *i ) { // Not a duplicate

      if ( required_unique_by && e != i ) {
        const unsigned attr_hi_to_lo = hi_to_lo.attribute();
        const unsigned attr_existing = i->attribute();

        if ( attr_hi_to_lo == attr_existing ) {
          throw_required_unique( e_hi, e_lo, *i->entity(), identifier,
                                 required_unique_by);
        }
      }

      e_hi.m_connect.insert( i , hi_to_lo );
    }
  }

  {
    const ConnectSet::iterator e = e_lo.m_connect.end();
          ConnectSet::iterator i = e_lo.m_connect.begin();

    const Connect lo_to_hi( e_hi , identifier );

    i = std::lower_bound( i , e , lo_to_hi , LessConnect() );

    if ( e == i || lo_to_hi != *i ) { // Not a duplicate
      e_lo.m_connect.insert( i , lo_to_hi );
    }
  }
}

void Mesh::declare_connection( Entity & e ,
                               const std::vector<Connect> & con ,
                               const char * required_unique_by )
{
  for ( std::vector<Connect>::const_iterator
        i = con.begin() ; i != con.end() ; ++i ) {

    const unsigned r_ident  =   i->identifier();
          Entity & r_entity = * i->entity();

    declare_connection( e , r_entity , r_ident , required_unique_by );
  }
}

void Mesh::destroy_connection( Entity & e1 , Entity & e2 )
{
  ConnectSet::iterator i ;

  for ( i = e1.m_connect.end() ; i != e1.m_connect.begin() ; ) {
    --i ;
    if ( & e2 == i->entity() ) { i = e1.m_connect.erase( i ); }
  }

  for ( i = e2.m_connect.end() ; i != e2.m_connect.begin() ; ) {
    --i ;
    if ( & e1 == i->entity() ) { i = e2.m_connect.erase( i ); }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

