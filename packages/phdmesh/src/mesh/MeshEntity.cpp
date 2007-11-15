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

const char * connect_type_name( ConnectType r )
{
  static const char * connect_names[] = { "USES" , "USEDBY" , "ANONYMOUS" };

  return connect_names[ (unsigned) r ];
}

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

Connect::Connect( Entity & e , ConnectType t , unsigned id )
  : m_attr( attribute( e.entity_type(), t, id) ) , m_entity(&e)
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
operator << ( std::ostream & os , const Connect & con )
{
  const char * const name = connect_type_name( con.type() );
  const unsigned     id   = con.identifier();

  os << name << "[" << id << "]->" ;

  if ( con.entity() != NULL ) { print_entity_key( os , con.entity()->key() ); }

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

ConnectSpan Entity::connections( EntityType et ) const
{
  const unsigned lo_key = Connect::attribute( et ,   ConnectType(0) , 0 );
  const unsigned hi_key = Connect::attribute( et+1 , ConnectType(0) , 0 );

  ConnectSet::const_iterator i = m_connect.begin();
  ConnectSet::const_iterator e = m_connect.end();

  i = std::lower_bound( i , e , lo_key , LessConnectAttr() );
  e = std::lower_bound( i , e , hi_key , LessConnectAttr() );

  return ConnectSpan( i , e );
}

ConnectSpan Entity::connections( EntityType et , ConnectType t ) const
{
  const unsigned lo_key = Connect::attribute( et , t ,   0 );
  const unsigned hi_key = Connect::attribute( et , t+1 , 0 );

  ConnectSet::const_iterator i = m_connect.begin();
  ConnectSet::const_iterator e = m_connect.end();

  i = std::lower_bound( i , e , lo_key , LessConnectAttr() );
  e = std::lower_bound( i , e , hi_key , LessConnectAttr() );

  return ConnectSpan( i , e );
}

//----------------------------------------------------------------------

namespace {

inline bool equal_attr( const Connect & lhs , const Connect & rhs )
{
  const unsigned lhs_attr = lhs.attribute();
  const unsigned rhs_attr = rhs.attribute();
  return lhs_attr == rhs_attr ;
}

}

void Mesh::declare_connection_anon( Entity & e1 , Entity & e2 ,
                                    const unsigned identifier )
{
  const Connect e1_to_e2( e2 , Anonymous , identifier );

  const ConnectSet::iterator e = e1.m_connect.end();
        ConnectSet::iterator i = e1.m_connect.begin();

  i = std::lower_bound( i , e , e1_to_e2 , LessConnect() );

  if ( e == i || e1_to_e2 != *i ) {
    e1.m_connect.insert( i , e1_to_e2 );
  }
}

void Mesh::declare_connection( Entity & e1 , Entity & e2 ,
                               const unsigned identifier ,
                               const char * required_unique_by )
{
  static const char method_name[] = "phdmesh::Mesh::declare_connection" ;

  std::ostringstream msg ;

  const EntityType e1_type = e1.entity_type();
  const EntityType e2_type = e2.entity_type();

  Entity & e_hi = e1_type == e2_type ? e1 : ( e1_type < e2_type ? e2 : e1 );
  Entity & e_lo = e1_type == e2_type ? e2 : ( e1_type < e2_type ? e1 : e2 );

  const Connect hi_to_lo( e_lo , Uses , identifier );
  const Connect lo_to_hi( e_hi , UsedBy , identifier );

  {
    const ConnectSet::iterator e = e_hi.m_connect.end();
          ConnectSet::iterator i = e_hi.m_connect.begin();

    i = std::lower_bound( i , e , hi_to_lo , LessConnect() );

    if ( required_unique_by && e != i && equal_attr( hi_to_lo , *i ) ) {
      msg << method_name << "( " ;
      print_entity_key( msg , e1.key() );
      msg << " , " ;
      print_entity_key( msg , e2.key() );
      msg << " , " ;
      msg << identifier ;
      msg << " , " ;
      msg << required_unique_by ;
      msg << " ) FAILED : ALREADY HAS " ;
      print_entity_key( msg , e_hi.key() );
      msg << "->{ Uses , " ;
      msg << identifier ;
      msg << " }->" ;
      print_entity_key( msg , i->entity()->key() );

      throw std::invalid_argument(msg.str());
    }

    if ( e == i || hi_to_lo != *i ) {
      e_hi.m_connect.insert( i , hi_to_lo );
    }
  }

  {
    const ConnectSet::iterator e = e_lo.m_connect.end();
          ConnectSet::iterator i = e_lo.m_connect.begin();

    i = std::lower_bound( i , e , lo_to_hi , LessConnect() );

    if ( e == i || lo_to_hi != *i ) {
      e_lo.m_connect.insert( i , lo_to_hi );
    }
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

