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
#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>

#include <util/ValueIOS.hpp>

namespace phdmesh {

namespace {

struct compare_type_info {

  bool operator()( const std::type_info & lhs ,
                   const std::type_info & rhs ) const
  { return lhs.before( rhs ) && lhs != rhs ; }

  bool operator()( const ValueIOS<void> * const lhs ,
                   const std::type_info &       rhs ) const
  { return operator()( lhs->type() , rhs ); }

  compare_type_info() {}
}; 

}

//----------------------------------------------------------------------

const ValueIOS<void> *
ValueIOSPolicy::get_void( const std::type_info & scalar_type ) const
{
  const compare_type_info compare ;

  std::vector< const ValueIOS<void> *>::const_iterator
    i = std::lower_bound( m_ios.begin() , m_ios.end() , scalar_type , compare );

  return i != m_ios.end() && (*i)->type() == scalar_type ? *i : NULL ;
}

void ValueIOSPolicy::replace( const ValueIOS<void> & ios )
{
  const compare_type_info compare ;

  std::vector< const ValueIOS<void> *>::iterator
    i = std::lower_bound( m_ios.begin() , m_ios.end() , ios.type() , compare );

  if ( i != m_ios.end() && (*i)->type() == ios.type() ) {
    *i = & ios ;
  }
  else {
    const ValueIOS<void> * const tmp = & ios ;
    m_ios.insert( i , tmp );
  }
}

//----------------------------------------------------------------------

ValueIOSPolicy::~ValueIOSPolicy() {}

ValueIOSPolicy::ValueIOSPolicy() : m_ios()
{
  static const ValueIOS<short>          ios_short ;
  static const ValueIOS<unsigned short> ios_unsigned_short ;
  static const ValueIOS<int>            ios_int ;
  static const ValueIOS<unsigned int>   ios_unsigned_int ;
  static const ValueIOS<long>           ios_long ;
  static const ValueIOS<unsigned long>  ios_unsigned_long ;
  static const ValueIOS<float>          ios_float ;
  static const ValueIOS<double>         ios_double ;
  static const ValueIOS<std::string>    ios_string ;

  replace( ios_short );
  replace( ios_unsigned_short );
  replace( ios_int );
  replace( ios_unsigned_int );
  replace( ios_long );
  replace( ios_unsigned_long );
  replace( ios_float );
  replace( ios_double );
  replace( ios_string );
}

ValueIOSPolicy::ValueIOSPolicy( const ValueIOSPolicy & p )
  : m_ios( p.m_ios )
{}

//----------------------------------------------------------------------

void ValueIOS<short>
  ::get( std::istream & s, short & v ) const
{ s >> v ; }

void ValueIOS<short int>
  ::put( std::ostream & s, unsigned, const short int & v) const
{ s << v ; }

void ValueIOS<short int>
  ::tell( std::ostream & s, unsigned, const short int & ) const
{ s << "short" ; }

void ValueIOS<short int>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<short int>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<short int>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<unsigned short>
  ::get( std::istream & s, unsigned short & v ) const
{ s >> v ; }

void ValueIOS<unsigned short>
  ::put( std::ostream & s, unsigned, const unsigned short & v) const
{ s << v ; }

void ValueIOS<unsigned short>
  ::tell( std::ostream & s, unsigned, const unsigned short &  ) const
{ s << "unsigned short" ; }

void ValueIOS<unsigned short>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned short>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned short>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<int>
  ::get( std::istream & s, int & v ) const
{ s >> v ; }

void ValueIOS<int>
  ::put( std::ostream & s, unsigned, const int & v) const
{ s << v ; }

void ValueIOS<int>
  ::tell( std::ostream & s, unsigned, const int &  ) const
{ s << "int" ; }

void ValueIOS<int>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<int>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<int>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<unsigned>
  ::get( std::istream & s, unsigned & v ) const
{ s >> v ; }

void ValueIOS<unsigned>
  ::put( std::ostream & s, unsigned, const unsigned & v) const
{ s << v ; }

void ValueIOS<unsigned>
  ::tell( std::ostream & s, unsigned, const unsigned &  ) const
{ s << "unsigned" ; }

void ValueIOS<unsigned>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<long>
  ::get( std::istream & s, long & v ) const
{ s >> v ; }

void ValueIOS<long>
  ::put( std::ostream & s, unsigned, const long & v) const
{ s << v ; }

void ValueIOS<long>
  ::tell( std::ostream & s, unsigned, const long &  ) const
{ s << "long" ; }

void ValueIOS<long>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<long>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<long>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<unsigned long>
  ::get( std::istream & s, unsigned long & v ) const
{ s >> v ; }

void ValueIOS<unsigned long>
  ::put( std::ostream & s, unsigned, const unsigned long & v) const
{ s << v ; }

void ValueIOS<unsigned long>
  ::tell( std::ostream & s, unsigned, const unsigned long &  ) const
{ s << "unsigned long" ; }

void ValueIOS<unsigned long>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned long>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned long>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<float>
  ::get( std::istream & s, float & v ) const
{ s >> v ; }

void ValueIOS<float>
  ::put( std::ostream & s, unsigned, const float & v) const
{ s << v ; }

void ValueIOS<float>
  ::tell( std::ostream & s, unsigned, const float &  ) const
{ s << "float" ; }

void ValueIOS<float>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<float>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<float>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS<double>
  ::get( std::istream & s, double & v ) const
{ s >> v ; }

void ValueIOS<double>
  ::put( std::ostream & s, unsigned, const double & v) const
{ s << v ; }

void ValueIOS<double>
  ::tell( std::ostream & s, unsigned, const double &  ) const
{ s << "double" ; }

void ValueIOS<double>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<double>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<double>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void enum_tell( std::ostream & os,
                const char * const indent ,
                const ValueIOS_Enum * en )
{
  os << "{" ;

  if ( ! en[1].m_name ) {
    os << " " ;
    os << en->m_name ;
    os << " = " ;
    os << en->m_value ;
    os << " " ;
  }
  else {
    while ( en->m_name ) {
      if ( indent ) { os << std::endl << indent ; }
      else          { os << " " ; }
      os << en->m_name ;
      os << " = " ;
      os << en->m_value ;
      os << " " ;
      ++en ;
      if ( ! en->m_name ) { os << "," ; }
    }
  }
  os << "}" ;
}

}

void enum_tell( std::ostream & os, unsigned indent,
                const char * name , const ValueIOS_Enum * en )
{
  static const char buf[] = "                                                                               " ;
  static const unsigned buf_len = sizeof(buf) - 1 ;

  const char * const b = buf + buf_len - indent ;

  os << "enum " << name << " " ;
  enum_tell( os , b , en );
}


long enum_value_of_name( const ValueIOS_Enum * en , const char * name )
{
  const ValueIOS_Enum * i = en ;
  while ( i->m_name && strcasecmp( i->m_name , name ) ) { ++i ; }

  if ( ! i->m_name ) {
    std::ostringstream msg ;
    msg << name ;
    msg << "' NOT IN " ;
    enum_tell( msg , NULL , en );
    throw std::runtime_error( msg.str() );
  }
  return i->m_value ;
}

const char * enum_name_of_value( const ValueIOS_Enum * en , const long v )
{
  const ValueIOS_Enum * i = en ;

  while ( i->m_name && v != i->m_value ) { ++i ; }

  if ( ! i->m_name ) {
    std::ostringstream msg ;
    msg << v ;
    msg << "' NOT IN " ;
    enum_tell( msg , NULL , en );
    throw std::runtime_error( msg.str() );
  }
  return i->m_name ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

const ValueIOS_Enum * named_bool_values()
{
  static const ValueIOS_Enum values[] = {
    { "false" , false } ,
    { "true" , true } ,
    { NULL , 0 }
  };
  return values ;
}

}

void ValueIOS<bool>
  ::get( std::istream & s , bool & v ) const
{
  std::string name ; s >> name ;
  long tmp = enum_value_of_name( named_bool_values() , name.c_str() );
  v = bool( tmp );
}

void ValueIOS<bool>
  ::put( std::ostream & s , unsigned, const bool & v ) const
{ s << enum_name_of_value( named_bool_values() , (long) v ); }

void ValueIOS<bool>
  ::tell( std::ostream & s , unsigned indent , const bool & ) const
{ enum_tell( s , indent , "bool" , named_bool_values() ); }

void ValueIOS<bool>
  ::getp( std::istream & s , void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<bool>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<bool>
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void ValueIOS< std::string >
  ::get( std::istream & s , std::string & v ) const
{ s >> v ; }

void ValueIOS< std::string >
  ::put( std::ostream & s , unsigned, const std::string & v ) const
{ s << v ; }

void ValueIOS< std::string >
  ::tell( std::ostream & s , unsigned, const std::string & ) const
{ s << "<non-white-space-string>" ; }

void ValueIOS< std::string >
  ::getp( std::istream & s , void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::string >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::string >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS_Quoted::get( std::istream & s , std::string & v ) const
{
  v.clear();
  while ( s.good() && isspace( s.peek() ) ) { s.get(); }
  if ( s.peek() == '"' ) {
    s.get();
    while ( s.good() && s.peek() != '"' ) {
      char tmp[2] = { (char) s.get() , 0 };
      v.append( tmp );
    }
    if ( s.get() != '"' ) {
      std::string msg("Reading unterminated quoted string" );
      throw std::runtime_error(msg);
    }
  }
}

void ValueIOS_Quoted
  ::put( std::ostream & s, unsigned, const std::string & v ) const
{ s << "\"" << v << "\"" ; }

void ValueIOS_Quoted
  ::tell( std::ostream & s, unsigned, const std::string & ) const
{ s << "\"<non-quote-char-string>\"" ; }

//----------------------------------------------------------------------

}

