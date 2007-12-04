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

#include <elem_top/Topology.hpp>

namespace phdmesh {
namespace element {

//----------------------------------------------------------------------

template<>
const Topology & topology< Point , 0 >()
{
  StaticAssert< Point::number_vertex == 0 &&
                Point::number_edge   == 0 &&
                Point::number_side   == 0 &&
                Point::boundary >::ok();
  static const Topology top( Point() , 0 , 0 );
  return top ;
}

template<>
const Topology & topology< Sphere , 1 >()
{
  StaticAssert< Sphere::number_vertex == 1 &&
                Sphere::number_edge   == 0 &&
                Sphere::number_side   == 0 &&
                ! Sphere::boundary >::ok();
  static const Topology top( Sphere() , 1 , 1 );
  return top ;
}

}
}

