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

#ifndef phdmesh_HexahedronTopology_hpp
#define phdmesh_HexahedronTopology_hpp

#include <element/LocalTopology.hpp>
#include <element/QuadrilateralTopology.hpp>

namespace phdmesh {

template< unsigned NumNode = 0 > struct Hexahedron ;         // 3D only

//----------------------------------------------------------------------
/*--------------------------------------------------------------------*/
/**
 *  Linear 8-Node Hexahedron Nodes
 *
 *         7                    6
 *          o------------------o
 *         /|                 /|
 *        / |                / |
 *       /  |               /  |
 *      /   |              /   |
 *     /    |             /    |
 *    /     |            /     |
 * 4 /      |         5 /      |
 *  o------------------o       |
 *  |       |          |       |
 *  |     3 o----------|-------o 2
 *  |      /           |      /
 *  |     /            |     /
 *  |    /             |    /
 *  |   /              |   /
 *  |  /               |  /
 *  | /                | /
 *  |/                 |/
 *  o------------------o
 * 0                    1
 *
 *--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/**
 *  Quadratic 20-Node Hexahedron Nodes
 *
 *          7         18         6
 *           o--------o---------o
 *          /|                 /|
 *         / |                / |
 *        /  |               /  |
 *     19o   |            17o   |
 *      /  15o             /    o14
 *     /     |            /     |
 *  4 /      | 16        /      |
 *   o---------o--------o 5     |
 *   |       |       10 |       |
 *   |     3 o-------o--|-------o 2
 *   |      /           |      /
 *   |     /            |     /
 * 12o    /             o13  /
 *   |   o11            |   o9
 *   |  /               |  /
 *   | /                | /
 *   |/                 |/
 *   o---------o--------o
 *  0          8         1
 *
 *--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/**
 *  Quadratic 27-Node Hexahedron Nodes
 *
 *           x--------x---------x
 *          /|                 /|
 *         / |                / |
 *        /  |   22          /  |
 *       x   |    o         x   |
 *      /    x       o26   /    x     Node #20 is at centroid of element
 *     /     |            /     |
 *    /      |           /      |     "2D surface" containing nodes
 *   x---------x--------x       |      0,1,5,4 has node 25 at center....
 *   | 23o   |          |   o24 |
 *   |       x-------x--|-------x
 *   |      /           |      /
 *   |     /  25        |     /
 *   x    /    o        x    /
 *   |   x        o21   |   x
 *   |  /               |  /
 *   | /                | /
 *   |/                 |/
 *   x---------x--------x
 *
 *--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef
  MakeTypeList< IndexList< 0 , 1 ,  8 , 27 , 28 > ,
                IndexList< 1 , 2 ,  9 , 29 , 30 > ,
                IndexList< 2 , 3 , 10 , 31 , 32 > ,
                IndexList< 3 , 0 , 11 , 33 , 34 > ,
                IndexList< 4 , 5 , 16 , 43 , 44 > ,
                IndexList< 5 , 6 , 17 , 45 , 46 > ,
                IndexList< 6 , 7 , 18 , 47 , 48 > ,
                IndexList< 7 , 4 , 19 , 49 , 50 > ,
                IndexList< 0 , 4 , 12 , 35 , 39 > ,
                IndexList< 1 , 5 , 13 , 36 , 40 > ,
                IndexList< 2 , 6 , 14 , 37 , 41 > ,
                IndexList< 3 , 7 , 15 , 38 , 42 > >::type
  HexahedronEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 5 , 4 ,  8, 13, 16, 12,  25 > ,
                IndexList< 1 , 2 , 6 , 5 ,  9, 14, 17, 13,  24 > ,
                IndexList< 2 , 3 , 7 , 6 , 10, 15, 18, 14,  26 > ,
                IndexList< 0 , 4 , 7 , 3 , 12, 19, 15, 11,  23 > ,
                IndexList< 0 , 3 , 2 , 1 , 11, 10,  9,  8,  21 > ,
                IndexList< 4 , 5 , 6 , 7 , 16, 17, 18, 19,  22 > >::type
  HexahedronSideNodeMap ;

template<>
struct Hexahedron<0>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                8 , 0 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      HexahedronEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<> , 
                    BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> >::type ,
      HexahedronSideNodeMap >
{
  typedef Hexahedron<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Hexahedron<8>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                8 , 8 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      HexahedronEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<4> , 
                    BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> >::type ,
      HexahedronSideNodeMap >
{
  typedef Hexahedron<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Hexahedron<20>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                8 , 20 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      HexahedronEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<8> , 
                    BoundaryQuadrilateral<8> ,
                    BoundaryQuadrilateral<8> ,
                    BoundaryQuadrilateral<8> ,
                    BoundaryQuadrilateral<8> ,
                    BoundaryQuadrilateral<8> >::type ,
      HexahedronSideNodeMap >
{
  typedef Hexahedron<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Hexahedron<27>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                8 , 27 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      HexahedronEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<9> , 
                    BoundaryQuadrilateral<9> ,
                    BoundaryQuadrilateral<9> ,
                    BoundaryQuadrilateral<9> ,
                    BoundaryQuadrilateral<9> ,
                    BoundaryQuadrilateral<9> >::type ,
      HexahedronSideNodeMap >
{
  typedef Hexahedron<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

