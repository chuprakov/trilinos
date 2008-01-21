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
#include <limits>
#include <stdexcept>
#include <string>
#include <iostream>

#include <util/ValueIO.hpp>

namespace phdmesh {

namespace {

int peek_non_space( std::istream & s )
{
  while ( s.good() && isspace( s.peek() ) ) { s.get(); }
  return s.peek();
}

template<typename T>
size_t read_array( std::istream & s , T * v , size_t n )
{
  size_t i = 0 ;
  while ( i < n && ( s >> v[i] ) ) { ++i ; }
  s.clear( s.rdstate() & ~std::ios::failbit );
  return i ;
}

template<typename T>
void write_array( std::ostream & s , T * v , size_t n )
{
  for ( size_t i = 0 ; i < n ; ++i ) {
    if ( i ) { s << " " ; }
    s << v[i] ;
  }
}

void tell_range( std::ostream & s , size_t nget , size_t nput )
{
  const size_t max = std::numeric_limits<size_t>::max();
  if      ( 0   == nput && 1 == nget ) { s << " const" ; }
  else if ( 0   == nput )              { s << " const[" << nget << "]" ; }
  else if ( max == nput )              { s << "[*]" ; }
  else if ( 1   <  nput )              { s << "[" << nput << "]" ; }
}

}

//----------------------------------------------------------------------

size_t ValueIO<signed char>
  ::read( std::istream & s, signed char * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<signed char>
  ::write( std::ostream & s, const signed char * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<signed char>
  ::tell( std::ostream & s, const signed char * , size_t nget , size_t nput )
{
  s << "char" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<unsigned char>::read(
  std::istream & s, unsigned char * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<unsigned char>::write(
  std::ostream & s, const unsigned char * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<unsigned char>::tell(
  std::ostream & s, const unsigned char * , size_t nget , size_t nput )
{
  s << "unsigned char" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<short>::read( std::istream & s, short * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<short>::write( std::ostream & s, const short * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<short>
  ::tell( std::ostream & s, const short * , size_t nget , size_t nput )
{
  s << "short" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<unsigned short>::read(
  std::istream & s, unsigned short * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<unsigned short>::write(
  std::ostream & s, const unsigned short * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<unsigned short>::tell(
  std::ostream & s, const unsigned short * , size_t nget , size_t nput )
{
  s << "unsigned short" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<int>::read( std::istream & s, int * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<int>::write( std::ostream & s, const int * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<int>
  ::tell( std::ostream & s, const int * , size_t nget , size_t nput )
{
  s << "int" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<unsigned int>
  ::read( std::istream & s, unsigned int * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<unsigned int>
  ::write( std::ostream & s, const unsigned int * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<unsigned int>
  ::tell( std::ostream & s, const unsigned int * , size_t nget , size_t nput )
{
  s << "unsigned int" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<long>::read( std::istream & s, long * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<long>::write( std::ostream & s, const long * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<long>
  ::tell( std::ostream & s, const long * , size_t nget , size_t nput )
{
  s << "long" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<unsigned long>
  ::read( std::istream & s, unsigned long * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<unsigned long>
  ::write( std::ostream & s, const unsigned long * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<unsigned long>
  ::tell( std::ostream & s, const unsigned long * , size_t nget , size_t nput )
{
  s << "unsigned long" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<float>::read( std::istream & s, float * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<float>::write( std::ostream & s, const float * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<float>
  ::tell( std::ostream & s, const float * , size_t nget , size_t nput )
{
  s << "float" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO<double>::read( std::istream & s, double * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO<double>::write( std::ostream & s, const double * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO<double>
  ::tell( std::ostream & s, const double * , size_t nget , size_t nput )
{
  s << "double" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

std::string read_name( std::istream & s )
{
  std::string name ;
  int c = peek_non_space(s);

  if ( isalpha( c ) ) {
    while ( isalnum( c ) || c == '_' ) {
      const char tmp[2] = { (char) c , 0 };
      name.append( tmp );
      s.get();
      c = s.peek();
    }
  }

  return name ;
}

const ValueIO_Enum * find(
const ValueIO_Enum * i , const char * name )
{
  while ( i->name && strcasecmp( i->name , name ) ) { ++i ; }

  if ( ! i->name ) { i = NULL ; }

  return i ;
}

const ValueIO_Enum * find(
const ValueIO_Enum * i , long value )
{
  while ( i->name && i->value != value ) { ++i ; }

  if ( ! i->name ) { i = NULL ; }

  return i ;
}

const ValueIO_Enum * read(
const ValueIO_Enum * i , std::istream & s )
{
  const std::string name = read_name( s );
  return find( i , name.c_str() );
}

void tell( const ValueIO_Enum * i , std::ostream & s )
{
  s << "{" ;
  while ( i->name ) {
    s << " " << i->name ;
    ++i ;
    if ( i->name ) { s << " ," ; }
  }
  s << " }" ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

const ValueIO_Enum * named_bool_values()
{
  static const ValueIO_Enum values[] = {
    { "false" , false } ,
    { "true" , true } ,
    { NULL , 0 }
  };
  return values ;
}

}

size_t ValueIO<bool>::read( std::istream & s , bool * v , size_t n )
{
  const ValueIO_Enum * j ;
  size_t i = 0 ;
  while ( i < n && NULL != ( j = phdmesh::read( named_bool_values() , s ) ) ) {
    v[i] = j->value ;
  }
  return i ;
}

void ValueIO<bool>::write( std::ostream & s , const bool * v , size_t n )
{
  for ( size_t i = 0 ; i < n ; ++i ) {
    if ( i ) { s << " " ; }
    const ValueIO_Enum * j = find( named_bool_values(), (long) v[i] );
    if ( j ) { s << j->name ; }
    else     { s << "ERROR" ; }
  }
}

void ValueIO<bool>
  ::tell( std::ostream & s , const bool * , size_t nget , size_t nput )
{
  s << "enum bool " ;
  phdmesh::tell( named_bool_values() , s );
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

size_t ValueIO< std::string >::read(
  std::istream & s , std::string * v , size_t n )
{ return read_array( s , v , n ); }

void ValueIO< std::string >::write(
  std::ostream & s , const std::string * v , size_t n )
{ write_array( s , v , n ); }

void ValueIO< std::string >
  ::tell( std::ostream & s , const std::string * , size_t nget , size_t nput )
{
  s << "std::string'noSpace" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

size_t ValueIO_Quoted<std::string>
  ::read( std::istream & s , std::string * v , size_t n )
{
  size_t i = 0 ;

  while ( i < n && peek_non_space( s ) == '"' ) {
    std::string & vs = v[i] ;
    vs.clear();

    s.get();
    while ( s.good() && s.peek() != '"' ) {
      char tmp[2] = { (char) s.get() , 0 };
      vs.append( tmp );
    }
    if ( s.get() != '"' ) {
      std::string msg("Reading unterminated quoted string" );
      throw std::runtime_error(msg);
    }
  }
  return i ;
}

void ValueIO_Quoted<std::string>
  ::write( std::ostream & s, const std::string * v , size_t n )
{
  for ( size_t i = 0 ; i < n ; ++i ) {
    if ( i ) { s << " " ; }
    s << "\"" << v[i] << "\"" ;
  }
}

void ValueIO_Quoted<std::string>
  ::tell( std::ostream & s, const std::string * , size_t nget , size_t nput )
{
  s << "\"std::string'noQuote\"" ;
  tell_range( s , nget , nput );
}

//----------------------------------------------------------------------

}

