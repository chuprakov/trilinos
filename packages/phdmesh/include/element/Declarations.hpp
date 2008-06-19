/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#ifndef phdmesh_element_Declarations_hpp
#define phdmesh_element_Declarations_hpp

#include <mesh/Types.hpp>
#include <element/LocalTopology.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Attach an element local topology to a Part.
 *  There is at most one element topology allowed.
 */
void set_part_local_topology( Part & , const LocalTopology * singleton );

/** Attach an element local topology to a Part.
 *  There is at most one element topology allowed.
 */
template< class ElementTraits >
void set_part_local_topology( Part & p )
{ return set_part_local_topology( p , ElementTraits::descriptor() ); }

const LocalTopology * get_part_local_topology( Part & );

//----------------------------------------------------------------------
/** Declare an element with nodes conformal to the given local topology. */
Entity & declare_element( Mesh & mesh ,
                          const LocalTopology & ,
                          const unsigned elem_id ,
                          const unsigned node_id[] );

/** Declare an element with nodes conformal to the given local topology. */
template< class ElementTraits >
Entity & declare_element( Mesh & mesh ,
                          const unsigned elem_id ,
                          const unsigned node_id[] )
{
  return declare_element( mesh , * ElementTraits::descriptor() ,
                          elem_id , node_id );
}

/** Declare an element member of a part with a local topology
 *  and nodes conformal to that topology.
 */
Entity & declare_element( Mesh & mesh ,
                          Part & part ,
                          const unsigned elem_id ,
                          const unsigned node_id[] );

//----------------------------------------------------------------------

}

#endif

