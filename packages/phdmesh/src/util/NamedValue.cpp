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

bool value_ios_get_token( std::istream & s , int t , bool required )
{
  static const char method_name[] = "phdmesh::value_ios_get_token" ;

  // while ( s.good() && isspace( s.peek() ) ) { s.get(); }
   while ( isspace( s.peek() ) ) { s.get(); }

  const bool found = s.peek() == t ;

  if ( found ) {
    s.get();
  }
  else if ( required ) {
    const char c[2] = { (char) t , 0 };
    std::string msg ;
    msg.append( method_name )
       .append( " FAILED to find '" )
       .append( c )
       .append( "'" );
    throw std::runtime_error(msg);
  }
  return found ;
}

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

  bool operator()( const NamedValue<> * lhs ,
                   const char * rhs ) const
    { return (*compare)( lhs->name.c_str() , rhs ) < 0 ; }
};

std::vector< NamedValue<>* >::const_iterator
lower_bound( const std::vector<NamedValue<>* > & v , const char * n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset(compare_nocase) );
}

std::vector< NamedValue<>* >::iterator
lower_bound( std::vector<NamedValue<>*> & v , const char * n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset(compare_nocase) );
}

//----------------------------------------------------------------------

void verify_name_policy( const char * method_name , const char * name )
{
  if ( ! good_c_name( name ) ) {
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
                     const NamedValue<> & p ,
                     std::string path )
{
  if ( p.reference.type == typeid(NamedValueSet) ) {
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

      const std::vector<NamedValue<>*> & vec = ps->get_all();

      for ( std::vector<NamedValue<>*>::const_iterator
            i = vec.begin() ; i != vec.end() ; ++i ) {
        verify_no_loop( method , vs , **i , path );
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

NamedValue<void>::NamedValue( const char * n , const Reference<void> & r )
  : name(n), reference(r), m_sets() {}

//----------------------------------------------------------------------

NamedValueSet::~NamedValueSet()
{ clear(); }

NamedValueSet::NamedValueSet() : m_values()
{}

NamedValueSet::NamedValueSet( const NamedValueSet & ps ) : m_values()
{ assign( ps.m_values ); }

void NamedValueSet::clear()
{
  NamedValueSet * myself = this ;
  while ( ! m_values.empty() ) {
    NamedValue<> & v = * m_values.back();
    remove_this( v.m_sets , myself );
    m_values.pop_back();
  }
}

void NamedValueSet::assign( const std::vector<NamedValue<>*> & v )
{
  static const char method_name[] = "phdmesh::NamedValueSet::assign" ;

  for ( std::vector<NamedValue<>*>::const_iterator
        i = v.begin() ; i != v.end() ; ++i ) {
    verify_name_policy( method_name , (*i)->name.c_str() );
  }

  clear();

  m_values = v ;

  NamedValueSet * myself = this ;
  for ( std::vector< NamedValue<> * >::iterator
        i = m_values.begin() ; i != m_values.end() ; ++i ) {
    (*i)->m_sets.push_back( myself );
  }
}

//----------------------------------------------------------------------

NamedValue<> *
NamedValueSet::m_find( const std::type_info & t ,
                       const std::string & n ,
                       const char sep ) const
{
  static const char method_name[] = "phdmesh::NamedValueSet::find" ;

  NamedValue<> * v = NULL ;

  const std::string::size_type len = n.size();
  const std::string::size_type p   = sep ? n.find( sep ) : len ;

  if ( len == p || std::string::npos == p ) {
    // Local:
    const std::vector<NamedValue<>*>::const_iterator
      i = lower_bound( m_values , n.c_str() );

    if ( m_values.end() != i && n == (*i)->name ) {
      verify_type( method_name , n , (*i)->reference.type , t );
      v = *i ;
    }
  }
  else {
    // Nested:
    std::string local_name  = n.substr( 0 , p ); // before sep
    std::string nested_name = n.substr( p + 1 ); // after  sep

    const std::vector<NamedValue<>*>::const_iterator
      i = lower_bound( m_values , local_name.c_str() );

    verify_type( method_name , local_name ,
                 (*i)->reference.type , typeid(NamedValueSet) );

    NamedValueSet & nested =
      static_cast< NamedValue<NamedValueSet> *>( *i )->value ;

    v = nested.m_find( t , nested_name , sep );
  }

  return v ;
}
//----------------------------------------------------------------------

NamedValue<> *
NamedValueSet::m_insert( NamedValue<> & v , bool replace )
{
  static const char method_name[] = "phdmesh::NamedValueSet::replace" ;

  NamedValue<> * m = NULL ;

  verify_name_policy( method_name , v.name.c_str() );

  verify_no_loop( method_name , this , v , std::string() );

  std::vector<NamedValue<>*>::iterator
    i = lower_bound( m_values , v.name.c_str() );

  if ( m_values.end() != i && v.name == (*i)->name ) {

    verify_type( method_name, v.name, (*i)->reference.type, v.reference.type );

    if ( replace ) {
      if ( & v != *i ) {
        NamedValueSet * myself = this ;
        v.m_sets.push_back( myself );
        remove_this( (*i)->m_sets , myself );

        m = *i ;
        *i = & v ;
      }
    }
    else {
      m = *i ;
    }
  }
  else {
    NamedValueSet * myself = this ;
    v.m_sets.push_back( myself );
    m = & v ;
    m_values.insert( i , m );
    if ( replace ) { m = NULL ; }
  }

  return m ;
}

void NamedValueSet::remove( NamedValue<> & p )
{
  std::vector<NamedValue<>*>::iterator
    i = lower_bound( m_values , p.name.c_str() );

  if ( m_values.end() != i && *i == & p ) {
    m_values.erase( i );
    remove_this( p.m_sets , this );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

size_t ValueIO<NamedValueSet>::read(
  std::istream & s , NamedValueSet * v , size_t n )
{
  size_t i = 0 ;

  for ( ; i < n && value_ios_get_token( s , '{' , false ) ; ++i ) {
    NamedValueSet & nv = v[i] ;

    std::string ex_msg ;
    std::string name ;
    NamedValue<> * m ;

    while ( ! ex_msg.size() &&
            ! value_ios_get_token( s , '}' , false ) ) {

      name = read_name( s );

      if ( name.empty() ) {
        ex_msg.append(": failed to read <name>");
      }
      else if ( ! value_ios_get_token( s , '=' , false ) ) {
        ex_msg.append(": failed to read '='");
      }
      else if ( NULL == ( m = nv.find<void>( name ) ) ) {
        ex_msg.append(": '");
        ex_msg.append(name);
        ex_msg.append("' is not a member");
      }
      else if ( 0 == m->reference.read(s) ) {
        ex_msg.append(": failed to read a value for '");
        ex_msg.append(name);
        ex_msg.append("'");
      }
      else if ( ! value_ios_get_token( s , ';' , false ) ) {
        ex_msg.append(": failed to read ';'");
      }
    }

    if ( ex_msg.size() ) {
      std::string msg( "operator<<(std::istream &, phdmesh::NamedValueSet &)" );
      msg.append( ex_msg );
      throw std::runtime_error( msg );
    }
  }

  return i ;
}

void ValueIO<NamedValueSet>::write(
  std::ostream & s , const NamedValueSet * v , size_t n )
{
  for ( size_t j = 0 ; j < n ; ++j ) {
    const NamedValueSet & nv = v[j] ;

    const std::vector< NamedValue<> * > & vec = nv.get_all();

    s << "{" ;

    if ( ! vec.empty() ) {
      s << std::endl ;

      for ( std::vector<NamedValue<>* >::const_iterator
            i = vec.begin() ; i != vec.end() ; ++i ) {

        s << (*i)->name ;
        s << " = " ;
        (*i)->reference.write( s );
        s << " ;" << std::endl ;
      }
    }

    s << "}" ;
  }
}

//----------------------------------------------------------------------

void ValueIO<NamedValueSet>::tell(
  std::ostream & s , const NamedValueSet * v , size_t nget , size_t )
{
  if ( v && nget ) {

    for ( size_t j = 0 ; j < nget ; ++j ) {

      const std::vector< NamedValue<>* > & vec = v[j].get_all();

      s << "{" ;

      if ( ! vec.empty() ) {

        s << std::endl ;

        for ( std::vector<NamedValue<>*>::const_iterator
              i = vec.begin() ; i != vec.end() ; ++i ) {
          NamedValue<> & nv = **i ;

          s << nv.name ;
          s << " = " ;
          nv.reference.tell( s );
          s << " ;" ;
          s << std::endl ;
        }
      }

      s << "}" ;
    }
  }
}

std::ostream & operator << ( std::ostream & s , const NamedValueSet & v )
{ ValueIO<NamedValueSet>::write( s , & v , 1 ); return s ; }

std::istream & operator >> ( std::istream & s , NamedValueSet & v )
{ ValueIO<NamedValueSet>::read( s , & v , 1 ); return s ; }

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


