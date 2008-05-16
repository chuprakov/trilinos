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

#include <util/Dimension.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

void stride_verify( unsigned n , const unsigned * stride , const char * who )
{
  bool ok = true ;
  for ( unsigned i = 0 ; ok && i < n ; ++i ) { ok = stride[i] ; }
  for ( unsigned i = 1 ; ok && i < n ; ++i ) {
    ok = ! ( stride[i] % stride[i-1] );
  }
  if ( ! ok ) {
    std::ostringstream msg ;
    msg << who ;
    msg << " HAS BAD DIMENSION STRIDE = {" ;
    for ( unsigned i = 0 ; i < n ; ++i ) { msg << stride[i] << " " ; }
    msg << "}" ;
    throw std::runtime_error( msg.str() );
  }
}
    
void stride_copy( unsigned n , unsigned * dst , const unsigned * src ,
                  const char * who )
{
  stride_verify( n , src , who );
  for ( unsigned i = 0 ; i < n ; ++i ) { dst[i] = src[i] ; }
}

unsigned stride_map( unsigned n , const unsigned * stride ,
                                  const unsigned * index )
{
  unsigned offset = n ? index[0] : 0 ;
  for ( unsigned i = 1 ; i < n ; ++i ) { offset += index[i] * stride[i-1] ; }
  return offset ;
}

void stride_inv( unsigned n , const unsigned * stride ,
                 unsigned offset , unsigned * index )
{
  for ( unsigned i = n - 1 ; i ; ) {
    index[i] = offset / stride[i-1] ;
    offset %= stride[i-1] ;
  }
  if ( n ) { index[0] = offset ; }
}

void stride_size( unsigned n , const unsigned * stride , unsigned * size )
{
  if ( n ) { size[0] = stride[0] ; }
  for ( unsigned i = 1 ; i < n ; ++i ) { size[i] = stride[i] / stride[i-1] ; }
}

unsigned stride_size( unsigned n , const unsigned * stride )
{ return n ? stride[n-1] : 1 ; }


//----------------------------------------------------------------------


const DimensionTraits * DimensionAny::descriptor()
{ static const DimensionAny self ; return & self ; }

const char * DimensionAny::name() const
{ static const char n[] = "DimensionAny" ; return n ; }

std::string DimensionAny::encode( unsigned size , unsigned index ) const
{
  std::ostringstream tmp ;

  if ( size <= index ) {
    tmp << "DimensionAny::encode( " << size << " , " << index << " ) ERROR" ;
    throw std::runtime_error( tmp.str() );
  }

  // Pad for alphabetic sorting.

  if ( 1000 < size && index < 1000 ) { tmp << '0' ; }
  if (  100 < size && index <  100 ) { tmp << '0' ; }
  if (   10 < size && index <   10 ) { tmp << '0' ; }
  tmp << index ;

  return tmp.str();
}

unsigned DimensionAny::decode( unsigned size ,
                               const std::string & label ) const
{
  const unsigned index = size ?  atoi( label.c_str() ) : 0 ;
  return index ;
}

} // namespace phdmesh

//----------------------------------------------------------------------


