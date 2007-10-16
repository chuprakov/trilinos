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

  while ( s.good() && isspace( s.peek() ) ) { s.get(); }

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
lower_bound( const std::vector<NamedValue<>* > & v ,
             const NamedValueSet::CompareFunction comp ,
             const char * n )
{
  return std::lower_bound( v.begin(), v.end(), n , less_pset(comp) );
}

std::vector< NamedValue<>* >::iterator
lower_bound( std::vector<NamedValue<>*> & v ,
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
                     const NamedValue<> & p ,
                     std::string path )
{
  if ( p.scalar_type() == typeid(NamedValueSet) ) {
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

void NamedValue<>::throw_type( const std::type_info & t ) const
{
  std::string msg ;
  msg.append( "phdmesh::NamedValue<>[scalar_typeid(" )
     .append( scalar_type().name() )
     .append( ") != argument_typeid(" )
     .append( t.name() )
     .append( ")" );
  throw std::runtime_error( msg );
}

//----------------------------------------------------------------------

NamedValue<>::~NamedValue()
{
  while ( ! m_sets.empty() ) {
    NamedValueSet & ps = *m_sets.back() ;
    ps.remove( *this );
  }
}

NamedValue<>::NamedValue( const char * n ) : name(n), m_sets() {}

//----------------------------------------------------------------------

NamedValueSet::~NamedValueSet()
{ clear(); }

NamedValueSet::NamedValueSet()
  : m_compare( compare_nocase ),
    m_good( good_c_name ),
    m_ios( ValueIOS<NamedValueSet>::default_policy() ),
    m_values()
{}

NamedValueSet::NamedValueSet( CompareFunction arg_compare ,
                              GoodFunction    arg_good ,
                              const ValueIOSPolicy * arg_ios )
  : m_compare( arg_compare ? arg_compare : compare_nocase ),
    m_good( arg_good ? arg_good : good_c_name ),
    m_ios( arg_ios ? * arg_ios : ValueIOS<NamedValueSet>::default_policy() ),
    m_values()
{}

NamedValueSet::NamedValueSet( const NamedValueSet & ps )
  : m_compare( ps.m_compare ),
    m_good(    ps.m_good ),
    m_ios(     ps.m_ios ),
    m_values()
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
    verify_name_policy( method_name , m_good , (*i)->name.c_str() );
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
      i = lower_bound( m_values , m_compare , n.c_str() );

    if ( m_values.end() != i && n == (*i)->name ) {
      verify_type( method_name , n , (*i)->scalar_type() , t );
      v = *i ;
    }
  }
  else {
    // Nested:
    std::string local_name  = n.substr( 0 , p ); // before sep
    std::string nested_name = n.substr( p + 1 ); // after  sep

    const std::vector<NamedValue<>*>::const_iterator
      i = lower_bound( m_values , m_compare , local_name.c_str() );

    verify_type( method_name , local_name ,
                 (*i)->scalar_type() , typeid(NamedValueSet) );

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

  verify_name_policy( method_name , m_good , v.name.c_str() );

  verify_no_loop( method_name , this , v , std::string() );

  std::vector<NamedValue<>*>::iterator
    i = lower_bound( m_values , m_compare , v.name.c_str() );

  if ( m_values.end() != i && v.name == (*i)->name ) {

    verify_type( method_name, v.name, (*i)->scalar_type(), v.scalar_type() );

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
    i = lower_bound( m_values , m_compare , p.name.c_str() );

  if ( m_values.end() != i && *i == & p ) {
    m_values.erase( i );
    remove_this( p.m_sets , this );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

ValueIOS<NamedValueSet>::ValueIOS() : ValueIOS<void>() {}

ValueIOS<NamedValueSet>::~ValueIOS() {}

const std::type_info & ValueIOS<NamedValueSet>::type() const
{ return typeid(NamedValueSet); }

void ValueIOS<NamedValueSet>::getp( std::istream & s , void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<NamedValueSet>::get( std::istream & s , NamedValueSet & v ) const
{
  static const char method[] = "phdmesh::ValueIOS<NamedValueSet>::get" ;

  if ( value_ios_get_token(s,'{',false) ) {
    std::string ex_msg ;
    std::string name ;

    try {
      NamedValue<> * m = NULL ;

      while ( ! value_ios_get_token( s , '}' , false ) ) {

        s >> name ;

        if ( ! s.good() ) {
          throw std::runtime_error( std::string("bad std::istream") );
        }

        value_ios_get_token( s , '=' , true );

        m = v.find<void>( name );

        if ( m == NULL ) {
          throw std::runtime_error( std::string("not found") );
        }

        const ValueIOS<void> * const ios = v.m_ios.get_void(m->scalar_type());

        if ( ios == NULL ) {
          throw std::runtime_error( std::string("no ValueIOS") );
        }

        for ( size_t i = 0 ; ! value_ios_get_token(s,';',false) ; ++i ) {
          void * const p = m->put_void(i);
          if ( ! p ) {
            throw std::runtime_error( std::string("too many values") );
          }
          ios->getp( s , p );
        }
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

  const char * const b  = buf + buf_len - indent ;
  const char * const b2 = buf + buf_len - ( indent + 2 );

  const std::vector< NamedValue<> * > & vec = v.get_all();

  if ( vec.empty() ) {
    s << "{}" ;
  }
  else {
    s << "{" ;

    for ( std::vector<NamedValue<>* >::const_iterator
          i = vec.begin() ; i != vec.end() ; ++i ) {

      s << std::endl << b2 << (*i)->name ;

      const ValueIOS<void> * const ios = v.m_ios.get_void((*i)->scalar_type());

      if ( ios ) {
        s << " = " ;
        const size_t get_max = (*i)->get_max();
        for ( unsigned j = 0 ; j < get_max ; ++j ) {
          ios->putp( s , indent + 2 , (*i)->get_void(j) );
          s << " " ;
        }
        s << ";" ;
      }
      else {
        s << " <?> ;" ;
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

  const char * const b  = buf + buf_len - indent ;
  const char * const b2 = buf + buf_len - ( indent + 2 );

  const std::vector< NamedValue<>* > & vec = v.get_all();

  if ( vec.empty() ) {
    s << "{}" ;
  }
  else {
    s << "{" ;

    for ( std::vector<NamedValue<>*>::const_iterator
          i = vec.begin() ; i != vec.end() ; ++i ) {

      s << std::endl << b2 << (*i)->name ;

      const ValueIOS<void> * const ios = v.m_ios.get_void((*i)->scalar_type());

      if ( ios ) {
        const size_t put_max = (*i)->put_max();
        if ( put_max == std::numeric_limits<size_t>::max() ) {
           s << "[*]" ;
        }
        else if ( 1 < put_max ) {
           s << "[" << put_max << "]" ;
        }
        s << " = " ;
        ios->tellp( s , indent + 2 , (*i)->get_void(0) );
        s << " ;" ;
      }
      else {
        s << " <?> ;" ;
      }
    }

    s << std::endl << b << "}" ;
  }
}

const ValueIOS<NamedValueSet> & ValueIOS<NamedValueSet>::singleton()
{ static const ValueIOS<NamedValueSet> io ; return io ; }

const ValueIOSPolicy & ValueIOS<NamedValueSet>::default_policy()
{
  static ValueIOSPolicy ios_policy ;
  static bool first = true ;

  if ( first ) { ios_policy.replace( singleton() ); }

  return ios_policy ;
}

std::ostream & operator << ( std::ostream & s , const NamedValueSet & v )
{
  const ValueIOS<NamedValueSet> * const io = v.m_ios.get<NamedValueSet>();
  io->put( s , 0 , v );
  return s ;
}

std::istream & operator >> ( std::istream & s , NamedValueSet & v )
{
  const ValueIOS<NamedValueSet> * const io = v.m_ios.get<NamedValueSet>();
  io->get( s , v );
  return s ;
}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------


