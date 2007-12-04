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

namespace {

const Topology::Orientation * edge5_orientation()
{
  StaticAssert< 2 == Edge::number_vertex >::ok();
  StaticAssert< 0 == Edge::number_edge >::ok();
  StaticAssert< 0 == Edge::number_side >::ok();
  StaticAssert<      Edge::boundary >::ok();
  StaticAssert<   Orientation<Edge,0>::polarity >::ok();
  StaticAssert< ! Orientation<Edge,1>::polarity >::ok();

  static const unsigned edge5_DP[] = { 0 , 1 , 2 , 3 , 4 };
  static const unsigned edge5_DM[] = { 1 , 0 , 2 , 4 , 3 };

  static const Topology::Orientation orientation[] =
    { { edge5_DP , 0 , true },
      { edge5_DM , 0 , false } };

  return orientation ;
}

}

template<>
const Topology & topology< Edge , 2 >()
{
  static const Topology top( Edge(), 2, 5, NULL, NULL, edge5_orientation() );
  return top ;
}

template<>
const Topology & topology< Edge , 3 >()
{
  static const Topology top( Edge(), 3, 5, NULL, NULL, edge5_orientation() );
  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< Line1D , 2 >()
{
  static const Topology top( Line1D(), 2, 5 );
  return top ;
}

template<>
const Topology & topology< Line1D , 3 >()
{
  static const Topology top( Line1D(), 3, 5 );
  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< LineBar , 2 >()
{
  static const Topology & edge_top = topology<Edge,2>();

  static const Topology::Boundary edge =
    { edge_top.orientation[0].node_map , & edge_top };

  static const Topology top( LineBar() , 2 , 5 , & edge );

  return top ;
}

template<>
const Topology & topology< LineBar , 3 >()
{
  static const Topology & edge_top = topology<Edge,3>();

  static const Topology::Boundary edge =
    { edge_top.orientation[0].node_map , & edge_top };

  static const Topology top( LineBar() , 3 , 5 , & edge );

  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< Line2DShell , 2 >()
{
  static const Topology & edge_top = topology<Edge,2>();

  static const Topology::Boundary edges[] =
    { { edge_top.orientation[0].node_map , & edge_top },
      { edge_top.orientation[1].node_map , & edge_top } };

  static const Topology top( LineBar() , 2 , 5 , edges , edges );

  return top ;
}

template<>
const Topology & topology< Line2DShell , 3 >()
{
  static const Topology & edge_top = topology<Edge,3>();

  static const Topology::Boundary edges[] =
    { { edge_top.orientation[0].node_map , & edge_top },
      { edge_top.orientation[1].node_map , & edge_top } };

  static const Topology top( LineBar() , 3 , 5 , edges , edges );

  return top ;
}

//----------------------------------------------------------------------

}
}

