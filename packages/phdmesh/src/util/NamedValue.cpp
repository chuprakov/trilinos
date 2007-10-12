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

#include <ctype.h>
#include <strings.h>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include <util/NamedValue.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

bool good_c_name( const char * name )
{
  bool result = isalpha( *name );

  while ( result && *++name ) {
    result = isalnum(*name) || *name == '_' ;
  }
  return result ;
}

int compare_nocase( const char * lhs , const char * rhs )
{
  return strcasecmp( lhs , rhs );
}

//----------------------------------------------------------------------

struct less_pset {

  int (*compare)( const char * , const char * );

  less_pset( int (*comp)( const char * , const char * ) ) : compare(comp) {}

  bool operator()( const NamedValueSet::MemberVoid & lhs ,
                   const char * rhs ) const
    { return (*compare)( lhs.first->name.c_str() , rhs ) < 0 ; }
};

std::vector< NamedValueSet::MemberVoid >::const_iterator
lower_bound( const std::vector<NamedValueSet::MemberVoid> & v ,
             const NamedValueSet::CompareFunction comp ,
             const char * n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset(comp) );
}

std::vector< NamedValueSet::MemberVoid >::iterator
lower_bound( std::vector<NamedValueSet::MemberVoid> & v ,
             const NamedValueSet::CompareFunction comp ,
             const char * n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset(comp) );
}

//----------------------------------------------------------------------

void verify_name_policy( const char * method_name ,
                         NamedValueSet::GoodFunction good ,
                         const char * name )
{
  if ( ! (*good)( name ) ) {
    std::string msg ;
    msg.append( method_name )
       .append( "( " )
       .append( name )
       .append( " ) FAILED: bad name" );
    throw std::runtime_error( msg );
  }
}

void verify_type( const char * method_name ,
                  const std::string    & name ,
                  const std::type_info & old_t ,
                  const std::type_info & new_t )
{
  if ( new_t != typeid(void) && new_t != old_t ) {
    std::string msg ;
    msg.append( method_name )
       .append( "( " )
       .append( name )
       .append( " ) FAILED: attempt to change type from typeid(" )
       .append( old_t.name() )
       .append( ") to typeid(" )
       .append( new_t.name() )
       .append( ")" );
    throw std::runtime_error( msg );
  }
}

void verify_no_loop( const char * const method ,
                     const NamedValueSet * const vs ,
                     const NamedValue<void> & p ,
                     std::string path )
{
  if ( p.value_type == typeid(NamedValueSet) ) {
    const NamedValueSet * ps =
      & static_cast< const NamedValue<NamedValueSet> & >(p).value ;

    if ( vs == ps ) {
      std::string msg ;
      msg.append( method )
         .append(" FAILED, Loop attempted: " )
         .append( p.name )
         .append( "." )
         .append( path )
         .append( p.name );
      throw std::runtime_error( msg );
    }
    else {
      path.append( p.name ).append( "." );

      const std::vector<NamedValueSet::MemberVoid> & vec = ps->get_all();

      for ( std::vector<NamedValueSet::MemberVoid>::const_iterator
            i = vec.begin() ; i != vec.end() ; ++i ) {
        verify_no_loop( method , vs , *(i->first) , path );
      }
    }
  }
}

//----------------------------------------------------------------------

void remove_this( std::vector< NamedValueSet *> & v , NamedValueSet * const ps )
{
  if ( ps ) {
    std::vector<NamedValueSet*>::iterator i ;
    for ( i = v.begin() ; i != v.end() && ps != *i ; ++i );
    if ( i != v.end() ) { v.erase( i ); }
  }
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

NamedValue<void>::~NamedValue()
{
  while ( ! m_sets.empty() ) {
    NamedValueSet & ps = *m_sets.back() ;
    ps.remove( *this );
  }
}

NamedValue<void>::NamedValue( const std::type_info & t , const char * n )
  : name(n), value_type(t), m_sets() {}

//----------------------------------------------------------------------

NamedValueSet::~NamedValueSet()
{
  for ( std::vector< MemberVoid >::iterator
        i = m_values.begin() ; i != m_values.end() ; ++i ) {
    remove_this( i->first->m_sets , this );
  }
}

NamedValueSet::NamedValueSet()
  : m_good( good_c_name ),
    m_compare( compare_nocase ),
    m_values()
{}

NamedValueSet::NamedValueSet( CompareFunction arg_compare ,
                              GoodFunction    arg_good )
  : m_good( arg_good ? arg_good : good_c_name ),
    m_compare( arg_compare ? arg_compare : compare_nocase ),
    m_values()
{}

NamedValueSet::NamedValueSet( const NamedValueSet & ps )
  : m_good(    ps.m_good ),
    m_compare( ps.m_compare ),
    m_values(  ps.m_values )
{
  NamedValueSet * myself = this ;
  for ( std::vector< MemberVoid >::iterator
        i = m_values.begin() ; i != m_values.end() ; ++i ) {
    i->first->m_sets.push_back( myself );
  }
}

void NamedValueSet::clear() { m_values.clear(); }

void NamedValueSet::assign( const std::vector<NamedValueSet::MemberVoid> & v )
{
  static const char method_name[] = "phdmesh::NamedValueSet::assign" ;

  for ( std::vector<MemberVoid>::const_iterator
        i = v.begin() ; i != v.end() ; ++i ) {
    verify_name_policy( method_name , m_good , i->first->name.c_str() );
  }
  m_values = v ;
}

//----------------------------------------------------------------------

NamedValueSet::MemberVoid
NamedValueSet::m_get( const std::type_info & t ,
                      const std::string & n ,
                      const char sep ) const
{
  static const char method_name[] = "phdmesh::NamedValueSet::get" ;

  NamedValue<void>     * vs = NULL ;
  const ValueIOS<void> * io = NULL ;

  MemberVoid m( vs , io );

  const std::string::size_type len = n.size();
  const std::string::size_type p   = sep ? n.find( sep ) : len ;

  if ( len == p || std::string::npos == p ) {
    // Local:
    const std::vector<MemberVoid>::const_iterator
      i = lower_bound( m_values , m_compare , n.c_str() );

    if ( m_values.end() != i && n == i->first->name ) {
      verify_type( method_name , n , i->first->value_type , t );
      m = *i ;
    }
  }
  else {
    // Nested:
    std::string local_name  = n.substr( 0 , p ); // before sep
    std::string nested_name = n.substr( p + 1 ); // after  sep

    const std::vector<MemberVoid>::const_iterator
      i = lower_bound( m_values , m_compare , local_name.c_str() );

    verify_type( method_name , local_name ,
                 i->first->value_type , typeid(NamedValueSet) );

    NamedValueSet & nested =
      static_cast< NamedValue<NamedValueSet> *>( i->first )->value ;

    m = nested.m_get( t , nested_name , sep );
  }

  return m ;
}
//----------------------------------------------------------------------

NamedValueSet::MemberVoid
NamedValueSet::m_insert( NamedValueSet::MemberVoid m )
{
  static const char method_name[] = "phdmesh::NamedValueSet::insert" ;

  verify_name_policy( method_name , m_good , m.first->name.c_str() );

  verify_no_loop( method_name , this , *m.first , std::string() );

  std::vector<MemberVoid>::iterator
    i = lower_bound( m_values , m_compare , m.first->name.c_str() );

  if ( m_values.end() != i && m.first->name == i->first->name ) {
    verify_type( method_name , m.first->name ,
                               i->first->value_type , m.first->value_type );
    m = *i ;
  }
  else {
    NamedValueSet * myself = this ;
    m.first->m_sets.push_back( myself );
    m_values.insert( i , m );
  }

  return m ;
}

NamedValueSet::MemberVoid
NamedValueSet::m_replace( NamedValueSet::MemberVoid m )
{
  static const char method_name[] = "phdmesh::NamedValueSet::replace" ;

  verify_name_policy( method_name , m_good , m.first->name.c_str() );

  verify_no_loop( method_name , this , *m.first , std::string() );

  std::vector<MemberVoid>::iterator
    i = lower_bound( m_values , m_compare , m.first->name.c_str() );

  if ( m_values.end() != i && m.first->name == i->first->name ) {

    verify_type( method_name , m.first->name ,
                               i->first->value_type , m.first->value_type );

    if ( m.first != i->first ) {
      NamedValueSet * myself = this ;
      m.first->m_sets.push_back( myself );
      remove_this( i->first->m_sets , myself );

      NamedValue<void> * tmp = i->first ;
      i->first = m.first ;
      m.first = tmp ;
    }
     
    if ( m.second ) {
      if ( m.second != i->second ) {
        const ValueIOS<void> * tmp = i->second ;
        i->second = m.second ;
        m.second = tmp ;
      }
    }
  }
  else {
    NamedValueSet * myself = this ;
    m.first->m_sets.push_back( myself );
    m_values.insert( i , m );
    m.first = NULL ;
    m.second = NULL ;
  }

  return m ;
}

void NamedValueSet::remove( NamedValue<void> & p )
{
  std::vector<MemberVoid>::iterator
    i = lower_bound( m_values , m_compare , p.name.c_str() );

  if ( m_values.end() != i && i->first == & p ) {
    m_values.erase( i );
    remove_this( p.m_sets , this );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

ValueIOS<NamedValueSet>::ValueIOS() : ValueIOS<void>() {}

ValueIOS<NamedValueSet>::~ValueIOS() {}

void ValueIOS<NamedValueSet>::getp( std::istream & s , void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<NamedValueSet>::get( std::istream & s , NamedValueSet & v ) const
{
  static const char method[] = "phdmesh::ValueIOS<NamedValueSet>::get" ;

  if ( value_ios_get_block_begin(s) ) {
    std::string ex_msg ;
    std::string name ;

    try {
      NamedValueSet::MemberVoid m ;

      while ( ! value_ios_get_block_end( s ) ) {

        s >> name ;

        if ( ! s.good() ) {
          throw std::runtime_error( std::string("bad std::istream") );
        }

        value_ios_get_token( s , '=' , true );

        m = v.get_ios<void>( name );

        if ( m.first == NULL ) {
          throw std::runtime_error( std::string("not found") );
        }

        if ( m.second == NULL ) {
          throw std::runtime_error( std::string("no ValueIOS") );
        }

        m.second->getp( s , m.first->pointer() );
      }
    }
    catch( const std::exception & x ) {
      ex_msg.append( x.what() );
    }
    catch( ... ) {
      ex_msg.append( "unknown exception" );
    }

    if ( ex_msg.size() ) {
      std::string msg ;
      msg.append( method )
         .append( " FAILED reading '" )
         .append( name )
         .append( "' due to:\n" )
         .append( ex_msg );
      throw std::runtime_error( msg );
    }
  }
}

void ValueIOS<NamedValueSet>::putp( std::ostream & s , unsigned indent ,
                                    const void * v ) const
{ put( s , indent, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<NamedValueSet>::put( std::ostream & s , unsigned indent ,
                                   const NamedValueSet & v ) const
{
  static const char buf[] = "                                                                               " ; 
  static const unsigned buf_len = sizeof(buf) - 1 ;

  const char * const b = buf + buf_len - indent ;

  const std::vector< NamedValueSet::MemberVoid > & vec = v.get_all();

  if ( vec.empty() ) {
    s << "{}" ;
  }
  else {
    s << "{" ;

    for ( std::vector<NamedValueSet::MemberVoid>::const_iterator
          i = vec.begin() ; i != vec.end() ; ++i ) {

      s << std::endl << b << i->first->name ;

      if ( i->second ) {
        s << " = " ;
        i->second->putp( s , indent + 2 , i->first->pointer() );
      }
      else {
        s << "< ValueIOS = NULL >" ;
      }
    }

    s << std::endl << b << "}" ;
  }
}

void ValueIOS<NamedValueSet>::tellp( std::ostream & s , unsigned indent ,
                                     const void * v ) const
{ tell( s , indent, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<NamedValueSet>::tell( std::ostream & s , unsigned indent ,
                                    const NamedValueSet & v ) const
{
  static const char buf[] = "                                                                               " ; 
  static const unsigned buf_len = sizeof(buf) - 1 ;

  const char * const b = buf + buf_len - indent ;

  const std::vector< NamedValueSet::MemberVoid > & vec = v.get_all();

  if ( vec.empty() ) {
    s << "{}" ;
  }
  else {
    s << "{" ;

    for ( std::vector<NamedValueSet::MemberVoid>::const_iterator
          i = vec.begin() ; i != vec.end() ; ++i ) {

      s << std::endl << b << i->first->name ;

      if ( i->second ) {
        s << " = " ;
        i->second->tellp( s , indent + 2 , i->first->pointer() );
      }
      else {
        s << "< ValueIOS = NULL >" ;
      }
    }

    s << std::endl << b << "}" ;
  }
}

const ValueIOS<NamedValueSet> & ValueIOS<NamedValueSet>::singleton()
{ static const ValueIOS<NamedValueSet> io ; return io ; }

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


