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
Connect::print( std::ostream & os ) const
{
  const char * const name = connect_type_name( type() );
  const unsigned     id   = identifier();

  os << name << "[" << id << "]->" ;

  if ( m_entity != NULL ) { print_entity_key( os , m_entity->key() ); }

  return os ;
}

std::ostream &
Entity::print( std::ostream & os , const std::string & lead ) const
{
  os << lead ;
  print_entity_key( os , key() );
  os << " Owner(P" << m_owner_rank << ") Connections {" ;
  
  const ConnectSet::const_iterator e = m_connect.end();
        ConnectSet::const_iterator i = m_connect.begin();

  for ( ; i != e ; ++i ) {
    os << std::endl << lead << "  " ;
    i->print( os );
  }
  os << " }" << std::endl ;

  return os ;
}

//----------------------------------------------------------------------

namespace {
const unsigned long zero_ul = 0 ;
}

Entity::Entity()
  : SetvMember<unsigned long>( zero_ul ),
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

  ConnectSpan result( m_connect.begin() , m_connect.end() );

  result.first  = std::lower_bound( result.first , result.second ,
                                    lo_key , LessConnectAttr() );

  result.second = std::lower_bound( result.first , result.second ,
                                    hi_key , LessConnectAttr() );

  return result ;
}

ConnectSpan Entity::connections( EntityType et , ConnectType t ) const
{
  const unsigned lo_key = Connect::attribute( et , t ,   0 );
  const unsigned hi_key = Connect::attribute( et , t+1 , 0 );

  ConnectSpan result( m_connect.begin() , m_connect.end() );

  result.first  = std::lower_bound( result.first , result.second ,
                                    lo_key , LessConnectAttr() );

  result.second = std::lower_bound( result.first , result.second ,
                                    hi_key , LessConnectAttr() );

  return result ;
}

//----------------------------------------------------------------------

void Entity::add_connection( const Connect & r ,
                             const char * required_unique_by )
{
  static const char method_name[] = "phdmesh::Entity::add_connection" ;

  std::ostringstream msg ;

  if ( r.entity() == NULL ) {
    msg << method_name << " ERROR : GIVEN NULL Entity" ;
    throw std::invalid_argument( msg.str() );
  }

  const ConnectSet::iterator e = m_connect.end();
        ConnectSet::iterator i = m_connect.begin();

  i = std::lower_bound( i , e , r , LessConnect() );

  bool flag = i != m_connect.end();

  if ( flag ) { // Not the end
    {
      const unsigned i_attr = i->attribute();
      const unsigned r_attr = r.attribute();

      flag = i_attr == r_attr ;
    }

    if ( flag ) { // Equal attributes

      {
        Entity * const i_entity = i->entity();
        Entity * const r_entity = r.entity();

        flag = i_entity == r_entity ;
      }

      // One or more connections with this attribute already exists

      if ( required_unique_by ) {
        msg << method_name << "( " ;
        r.print( msg );
        msg << " , " << required_unique_by
            << " ) ALREADY HAS " ;
        i->print( msg );

        throw std::invalid_argument(msg.str());
      }
    }
  }

  if ( ! flag ) {
    m_connect.insert( i , r );
  }
}

void Entity::remove_connections( Entity * const entity )
{
  ConnectSet::iterator i = m_connect.end();
  while ( i != m_connect.begin() ) {
    --i ;
    if ( entity == i->entity() ) {
      i = m_connect.erase( i );
    }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

