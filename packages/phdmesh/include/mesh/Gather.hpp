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

#include <util/Basics.hpp>
#include <mesh/Field.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

/** Gather data into an array T[ NValue X NConnect ]
 *  from connected entities.
 *
 *    double gdata[ 3 * NWork ];
 *    gather<3,NWork>( gdata , field , begin_entity , end_entity );
 */

template<unsigned NValue, unsigned NConnect>
class gather {
public:

  template<typename T,unsigned NDim>
  gather( T * dst ,
          const Field<T,NDim> & field ,
          const Entity ** i ,
          const Entity ** const j )
  {
    enum { NBlock = NValue * NConnect };

    for ( ; i < j ; ++i , dst += NBlock ) {
      ConnectSpan con = (*i)->connections( field.entity_type() );
      while ( con.first != con.second && con.identifier() < NConnect ) {
        const T * const s = con.first->entity()->data(field);
        if ( s ) { Copy<NValue>( dst + NValue * con.identifier() , s ); }
        ++con.first ;
      }
    }
  }
}

}

