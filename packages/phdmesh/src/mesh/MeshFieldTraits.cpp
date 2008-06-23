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

#include <stdlib.h>
#include <strings.h>

#include <sstream>
#include <stdexcept>
#include <mesh/FieldTraits.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const DimensionTag * Cartesian::descriptor()
{ static const Cartesian self ; return & self ; }

const char * Cartesian::name() const
{ static const char n[] = "Cartesian" ; return n ; }

std::string Cartesian::to_string( unsigned size , unsigned index ) const
{
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * label[] = { x , y , z };

  if ( size < 2 || 3 < size || size <= index ) {
    std::ostringstream msg ;
    msg << Cartesian::descriptor()->name();
    msg << " ERROR Size = " << size ;
    msg << " Index = " << index ;
    throw std::runtime_error( msg.str() );
  }

  return std::string( label[index] );
}

unsigned Cartesian::to_index( unsigned size , const std::string & arg ) const
{
  static const char x[] = "x" ;
  static const char y[] = "y" ;
  static const char z[] = "z" ;
  static const char * label[] = { x , y , z };

  unsigned index = size ;

  if ( 1 < size && size < 4 ) {
    const char * const c = arg.c_str();
    for ( index = 0 ; index < size && strcasecmp(c,label[index]) ; ++index );
  }

  if ( index == size ) {
    std::ostringstream msg ;
    msg << Cartesian::descriptor()->name();
    msg << " ERROR size = " << size ;
    msg << " label = " << arg ;
    throw std::runtime_error( msg.str() );
  }
  return index ;
}

//----------------------------------------------------------------------

const DimensionTag * Cylindrical::descriptor()
{ static const Cylindrical self ; return & self ; }

const char * Cylindrical::name() const
{ static const char n[] = "Cylindrical" ; return n ; }

std::string Cylindrical::to_string( unsigned size , unsigned index ) const
{
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * label[] = { r , a , z };

  if ( 3 < size || size <= index ) {
    std::ostringstream msg ;
    msg << Cylindrical::descriptor()->name();
    msg << " ERROR Size = " << size ;
    msg << " Index = " << index ;
    throw std::runtime_error( msg.str() );
  }

  return std::string( label[index] );
}

unsigned Cylindrical::to_index( unsigned size , const std::string & arg ) const
{
  static const char r[] = "r" ;
  static const char a[] = "a" ;
  static const char z[] = "z" ;
  static const char * label[] = { r , a , z };

  unsigned index = size ;

  if ( 1 < size && size < 4 ) {
    const char * const c = arg.c_str();
    for ( index = 0 ; index < size && strcasecmp(c,label[index]) ; ++index );
  }

  if ( index == size ) {
    std::ostringstream msg ;
    msg << Cylindrical::descriptor()->name();
    msg << " ERROR size = " << size ;
    msg << " label = " << arg ;
    throw std::runtime_error( msg.str() );
  }
  return index ;
}

//----------------------------------------------------------------------

}

