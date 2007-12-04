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
#include <stdlib.h>
#include <string.h>

#include <stdexcept>
#include <iostream>
#include <sstream>

#include <mesh/Field.hpp>
#include <mesh/Schema.hpp>
#include <mesh_io/FieldName.hpp>

namespace phdmesh {

namespace {

void internal_encode( std::ostringstream & arg_os ,
                      const unsigned       arg_ndim  ,
                      const unsigned * const arg_dim ,
                      const unsigned * const arg_indices )
{
  for ( unsigned i = arg_ndim ; i ; ) {
    --i ;
    const unsigned dim = arg_dim[i] ;
    const unsigned ind = arg_indices[i] ;

    // Pad for alphabetic sorting.

    arg_os << '_' ;
    if ( 1000 < dim && ind < 1000 ) { arg_os << '0' ; }
    if (  100 < dim && ind <  100 ) { arg_os << '0' ; }
    if (   10 < dim && ind <   10 ) { arg_os << '0' ; }
    arg_os << ind ;
  }
}

// Work right-to-left looking for digits separated by '_'

const char * internal_decode(
  const unsigned n ,
  unsigned * const arg_indices ,
  const char * const b , const char * c )
{
  bool result = true ;

  for ( unsigned i = 0 ; result && i < n ; ++i ) {

    result = false ;

    while ( ( c != b ) && ( result = isdigit(*--c) ) );

    if ( result && ( result = *c == '_' ) ) {
      ++c ;
      arg_indices[i] = atoi( c );
      --c ;
    }
  }
  return result ? c : NULL ;
}

}

FieldName::~FieldName() {}

FieldName::FieldName() {}

const char * FieldName::name() const
{
  static const char n[] = "phdmesh::FieldName" ;
  return n ;
}

void FieldName::verify( const unsigned arg_ndim  ,
                        const unsigned * const arg_dim ,
                        const unsigned * const arg_indices ) const
{
  std::ostringstream msg ;

  unsigned i ;

  for ( i = 0 ; i < arg_ndim && arg_indices[i] < arg_dim[i] ; ++i );

  if ( i < arg_ndim ) {
    msg << this->name() << " FAILED: Bad indices: {" ;
    for ( i = 0 ; i < arg_ndim ; ++i ) {
      msg << " " << arg_dim[i] ;
    }
    msg << " } <= {" ;
    for ( i = 0 ; i < arg_ndim ; ++i ) {
      msg << " " << arg_indices[i] ;
    }
    msg << " }" ;
    throw std::runtime_error( msg.str() );
  }
}

std::string FieldName::encode(
  const std::string &    arg_name ,
  const unsigned         arg_ndim ,
  const unsigned * const arg_dim ,
  const unsigned * const arg_indices ) const
{
  verify( arg_ndim , arg_dim , arg_indices );

  std::ostringstream tmp ;

  tmp << arg_name ;

  internal_encode( tmp , arg_ndim , arg_dim , arg_indices );

  return tmp.str();
}

bool FieldName::decode(
  const std::string & arg_text ,
        std::string & arg_name ,
  const unsigned      arg_ndim ,
        unsigned * const arg_indices ) const
{

  const char * const b = arg_text.c_str();
  const char * c = b + arg_text.size();

  bool result = NULL != ( c = internal_decode( arg_ndim, arg_indices, b, c ) );

  if ( result ) {
    const unsigned n = c - b ;
    arg_name.assign( b , n );
  }

  return result ;
}

//----------------------------------------------------------------------

namespace {

class FieldName_Tensor : public FieldName {
public:

  virtual std::string encode( const std::string & ,
                              const unsigned ,
                              const unsigned * const ,
                              const unsigned * const ) const ;

  virtual bool decode( const std::string & ,
                             std::string & ,
                       const unsigned ,
                             unsigned * const ) const ;

  virtual const char * name() const ;

  virtual ~FieldName_Tensor();

  FieldName_Tensor( const char * , unsigned , const char * const * );

private:
  const unsigned             m_num ;
  const char * const         m_name ;
  const char * const * const m_label ;

  FieldName_Tensor();
  FieldName_Tensor( const FieldName_Tensor & );
  FieldName_Tensor & operator = ( const FieldName_Tensor & );
};

FieldName_Tensor::~FieldName_Tensor() {}

FieldName_Tensor::FieldName_Tensor(
  const char * arg_name ,
  unsigned arg_num ,
  const char * const * arg_comp )
: m_num( arg_num ), m_name( arg_name ), m_label( arg_comp )
{ }

const char * FieldName_Tensor::name() const
{ return m_name ; }

std::string FieldName_Tensor::encode(
  const std::string &    arg_name ,
  const unsigned         arg_ndim ,
  const unsigned * const arg_dimension ,
  const unsigned * const arg_indices ) const
{
  verify( arg_ndim , arg_dimension , arg_indices );

  std::ostringstream tmp ;

  tmp << arg_name ;

  if ( 1 < arg_ndim ) {
    internal_encode( tmp , arg_ndim - 1, arg_dimension + 1, arg_indices + 1);
  }

  if ( arg_ndim ) {
    tmp << '_' << m_label[ arg_indices[0] ];
  }

  return tmp.str();
}

bool FieldName_Tensor::decode(
  const std::string & arg_text ,
        std::string & arg_name ,
  const unsigned      arg_ndim ,
        unsigned * const arg_indices ) const
{
  bool result = true ;

  // Work right-to-left looking for '_'

  const char * const b = arg_text.c_str();
  const char * c = b + arg_text.size();

  if ( arg_ndim ) {
    while ( b != c && *--c != '_' );
    result = *c == '_' ;

    if ( result ) {
      ++c ;
      unsigned index ;
      for ( index = 0 ; index < m_num && strcmp(c,m_label[index]) ; ++index );
      result = index < m_num ;
      arg_indices[0] = index ;
      --c ;
    }

    if ( 1 < arg_ndim ) {
      result = NULL != ( c = internal_decode( arg_ndim, arg_indices, b, c ) );
    }
  }

  if ( result ) {
    const unsigned n = c - b ;
    arg_name.assign( b , n );
  }

  return result ;
}

}

//----------------------------------------------------------------------

const FieldName & io_array_name()
{
  static const FieldName descriptor ;
  return descriptor ;
}


const FieldName & io_cartesian_vector()
{
  static const char name[] = "phdmesh::FieldName_Vector" ;
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * comp[] = { x , y , z };

  static const FieldName_Tensor descriptor( name , 3 , comp );

  return descriptor ;
}

const FieldName & io_cylindrical_vector()
{
  static const char name[] = "phdmesh::FieldName_Cylindrical" ;
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * comp[] = { r , a , z };

  static const FieldName_Tensor descriptor( name , 3 , comp );

  return descriptor ;
}

const FieldName * io_get_field_name( const Field<void,0> & f )
{
  Span< CSet::iterator<const FieldName> > attr =
    f.attributes().get<const FieldName>();

  return attr.empty() ? (const FieldName *) NULL : & *attr ;
}

void io_declare_field_name( Field<void,0> & f , const FieldName & d )
{
  static const char method[] = "phdmesh::io_declare" ;

  const FieldName * const exist = io_get_field_name( f );

  if ( ! exist ) {
    f.schema().declare_field_attribute<const FieldName>( f, & d , false );
  }
  else if ( exist != & d ) {
    std::ostringstream msg ;
    msg << method << "( " << f.name() << " , " << d.name()
        << " ) FAILED: Redeclaration from "
        << exist->name();
    throw std::logic_error( msg.str() );
  }
}

}


