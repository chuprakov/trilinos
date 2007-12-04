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

#ifndef phdmesh_ElementTopology_hpp
#define phdmesh_ElementTopology_hpp

#include <cstddef>
#include <elem_top/Shape.hpp>

//----------------------------------------------------------------------

namespace phdmesh {
namespace element {

class Topology ;

//----------------------------------------------------------------------

template< class Shape , unsigned NumNodes > const Topology & topology();

template<> const Topology & topology< Point , 0 >();
template<> const Topology & topology< Edge , 2 >();
template<> const Topology & topology< Edge , 3 >();
template<> const Topology & topology< TriangleFace , 3 >();
template<> const Topology & topology< TriangleFace , 6 >();
template<> const Topology & topology< QuadrilateralFace , 4 >();
template<> const Topology & topology< QuadrilateralFace , 8 >();
template<> const Topology & topology< QuadrilateralFace , 9 >();

template<> const Topology & topology< Sphere , 1 >();

template<> const Topology & topology< Line1D , 2 >();
template<> const Topology & topology< Line1D , 3 >();
template<> const Topology & topology< LineBar , 2 >();
template<> const Topology & topology< LineBar , 3 >();
template<> const Topology & topology< Line2DShell , 2 >();
template<> const Topology & topology< Line2DShell , 3 >();

template<> const Topology & topology< Triangle2D , 3 >();
template<> const Topology & topology< Triangle2D , 6 >();
template<> const Topology & topology< Quadrilateral2D , 4 >();
template<> const Topology & topology< Quadrilateral2D , 8 >();
template<> const Topology & topology< Quadrilateral2D , 9 >();

template<> const Topology & topology< TriangleShell , 3 >();
template<> const Topology & topology< TriangleShell , 6 >();
template<> const Topology & topology< QuadrilateralShell , 4 >();
template<> const Topology & topology< QuadrilateralShell , 8 >();
template<> const Topology & topology< QuadrilateralShell , 9 >();

template<> const Topology & topology< Tetrahedron , 4 >();
template<> const Topology & topology< Tetrahedron , 10 >();

template<> const Topology & topology< Pyramid , 5 >();
template<> const Topology & topology< Wedge , 6 >();

template<> const Topology & topology< Hexahedron , 8 >();
template<> const Topology & topology< Hexahedron , 20 >();
template<> const Topology & topology< Hexahedron , 27 >();

//----------------------------------------------------------------------
/** @class Topology
 *  Overlay nodes (a.k.a. points) on an element shape.
 *  Every vertex must have a node, these nodes must be ordered first and
 *  their ordering must be consistent with the shape's vertex ordering.
 */
class Topology {
public:

  /** Boundary edge and side members */

  struct Boundary {
    const unsigned * node_map ;
    const Topology * topology ;
  };

  /** Orientations of boundary shapes */

  struct Orientation {
    const unsigned * node_map ;
    unsigned short   rotation ;
    unsigned short   polarity ;
  };

  const unsigned number_vertex ;
  const unsigned number_edge ;
  const unsigned number_side ;
  const bool     boundary ;
  const unsigned number_orientation ;
  const unsigned number_node ;  // required nodes
  const unsigned maximum_node ; // understood nodes

  const Boundary    * const edge ;
  const Boundary    * const side ;
  const Orientation * const orientation ;

  template< unsigned NV , class EL , class SL ,
            unsigned MIN, unsigned MAX, bool B>
  Topology( const Shape<NV,EL,SL,MIN,MAX,B> ,
            unsigned            arg_number_node ,
            unsigned            arg_maximum_node ,
            const Boundary    * arg_edge = NULL ,
            const Boundary    * arg_side = NULL ,
            const Orientation * arg_orient = NULL )
    : number_vertex( Shape<NV,EL,SL,MIN,MAX,B>::number_vertex ),
      number_edge(   Shape<NV,EL,SL,MIN,MAX,B>::number_edge ),
      number_side(   Shape<NV,EL,SL,MIN,MAX,B>::number_side ),
      boundary(      Shape<NV,EL,SL,MIN,MAX,B>::boundary ),
      number_orientation(
        element::Orientation< Shape<NV,EL,SL,MIN,MAX,B> >::number_orientation ),
      number_node(  arg_number_node ),
      maximum_node( arg_maximum_node ),
      edge(         arg_edge ),
      side(         arg_side ),
      orientation(  arg_orient ) {}

  ~Topology() {}

  Topology( const Topology & rhs )
    : number_vertex(      rhs.number_vertex ),
      number_edge(        rhs.number_edge ),
      number_side(        rhs.number_side ),
      boundary(           rhs.boundary ),
      number_orientation( rhs.number_orientation ),
      number_node(        rhs.number_node ),
      maximum_node(       rhs.maximum_node ),
      edge(               rhs.edge ),
      side(               rhs.side ),
      orientation(        rhs.orientation ) {}

private:
  Topology();
  Topology & operator = ( const Topology & rhs );
};


} // namespace element
} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

