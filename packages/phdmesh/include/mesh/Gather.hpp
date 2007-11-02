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

#ifndef phdmesh_Gather_hpp
#define phdmesh_Gather_hpp

#include <iterator>

#include <util/Basics.hpp>
#include <mesh/Field.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

template<unsigned NValue, unsigned NConnect, typename T, unsigned NDim>
void gather_direct( T * dst ,
                    const Field<T,NDim> & field ,
                    const Entity ** i ,
                    const Entity ** const j )
{
  for ( ; i < j ; ++i ) {
    const ConnectSpan con = (*i)->connections( field.entity_type() );
    for ( unsigned k = 0 ; k < NConnect ; ++k , dst += NValue ) {
      const T * const s = con.first[k].entity()->data(field);
      Copy<NValue>( dst , s );
    }
  }
}

template<unsigned NValue, unsigned NConnect, typename T, unsigned NDim>
const Entity ** gather_direct_verify( const Field<T,NDim> & field ,
                                      const Entity ** i ,
                                      const Entity ** const j )
{
  enum { value_size = NValue * sizeof(T) };

  for ( bool result = true ; result && i < j ; ++i ) {
    ConnectSpan con = (*i)->connections( field.entity_type() );

    result = std::distance( con.first , con.second ) < NConnect ;

    for ( unsigned k = 0 ; result && k < NConnect ; ++k , ++con.first ) {
      result = k == con.first->identifier();

      if ( result ) {
        Entity & e = * con.first->entity();

        result = e.data(field);

        if ( result ) { result = e.data_size(field) == value_size ; }
      }
    }
  }
  return i ;
}

//----------------------------------------------------------------------

template<unsigned NValue, unsigned NConnect, typename T>
void gather_pointer( T * dst ,
                     const Field<T*,1> & field ,
                     const Entity * const * i ,
                     const Entity * const * const j )
{
  for ( ; i < j ; ++i ) {
    const T * const * const src = (*i)->data( field );
    for ( unsigned k = 0 ; k < NConnect ; ++k , dst += NValue ) {
      Copy<NValue>( dst , src[k] );
    }
  }
}

template<unsigned NConnect,typename T, unsigned N>
const Entity ** gather_pointer_setup( const Field<T,N>  & src_field ,
                                      const Field<T*,1> & ptr_field ,
                                      const Entity * const * i ,
                                      const Entity * const * const j )
{
  enum { value_size = NValue * sizeof(T) };
  enum { ptr_size   = NConnect * sizeof(T*) };

  for ( bool result = true ; result && i < j ; ++i ) {

    ConnectSpan con = (*i)->connections( field.entity_type() );

    result = std::distance( con.first , con.second ) < NConnect &&
             (*i)->data_size( ptr_field ) == ptr_size ;

    const T * const * const p = (*i)->data( ptr_field );

    for ( unsigned k = 0 ; result && k < NConnect ; ++k , ++con.first ) {
      result = k == con.first->identifier();

      if ( result ) {
        Entity & e = * con.first->entity();

        result = ( p[k] = e.data(field) );

        if ( result ) { result = e.data_size(field) == value_size ; }
      }
    }
  }
  return i ;
}

//----------------------------------------------------------------------

}

