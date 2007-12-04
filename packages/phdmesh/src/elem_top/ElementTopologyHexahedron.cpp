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

typedef Hexahedron Hex ;

namespace {

const unsigned * const * hex125_edge_map()
{
  StaticAssert< 0 == Hex::edge<0>::vertex<0>::ordinal >::ok();
  StaticAssert< 1 == Hex::edge<0>::vertex<1>::ordinal >::ok();
  StaticAssert< 1 == Hex::edge<1>::vertex<0>::ordinal >::ok();
  StaticAssert< 2 == Hex::edge<1>::vertex<1>::ordinal >::ok();
  StaticAssert< 2 == Hex::edge<2>::vertex<0>::ordinal >::ok();
  StaticAssert< 3 == Hex::edge<2>::vertex<1>::ordinal >::ok();
  StaticAssert< 3 == Hex::edge<3>::vertex<0>::ordinal >::ok();
  StaticAssert< 0 == Hex::edge<3>::vertex<1>::ordinal >::ok();
  StaticAssert< 4 == Hex::edge<4>::vertex<0>::ordinal >::ok();
  StaticAssert< 5 == Hex::edge<4>::vertex<1>::ordinal >::ok();
  StaticAssert< 5 == Hex::edge<5>::vertex<0>::ordinal >::ok();
  StaticAssert< 6 == Hex::edge<5>::vertex<1>::ordinal >::ok();
  StaticAssert< 6 == Hex::edge<6>::vertex<0>::ordinal >::ok();
  StaticAssert< 7 == Hex::edge<6>::vertex<1>::ordinal >::ok();
  StaticAssert< 7 == Hex::edge<7>::vertex<0>::ordinal >::ok();
  StaticAssert< 4 == Hex::edge<7>::vertex<1>::ordinal >::ok();
  StaticAssert< 0 == Hex::edge<8>::vertex<0>::ordinal >::ok();
  StaticAssert< 4 == Hex::edge<8>::vertex<1>::ordinal >::ok();
  StaticAssert< 1 == Hex::edge<9>::vertex<0>::ordinal >::ok();
  StaticAssert< 5 == Hex::edge<9>::vertex<1>::ordinal >::ok();
  StaticAssert< 2 == Hex::edge<10>::vertex<0>::ordinal >::ok();
  StaticAssert< 6 == Hex::edge<10>::vertex<1>::ordinal >::ok();
  StaticAssert< 3 == Hex::edge<11>::vertex<0>::ordinal >::ok();
  StaticAssert< 7 == Hex::edge<11>::vertex<1>::ordinal >::ok();

  static const unsigned en0[] = { 0,1,  8, 27,28 };
  static const unsigned en1[] = { 1,2,  9, 29,30 };
  static const unsigned en2[] = { 2,3, 10, 31,32 };
  static const unsigned en3[] = { 3,0, 11, 33,34 };
  static const unsigned en4[] = { 4,5, 16, 43,44 };
  static const unsigned en5[] = { 5,6, 17, 45,46 };
  static const unsigned en6[] = { 6,7, 18, 47,48 };
  static const unsigned en7[] = { 7,4, 19, 49,50 };
  static const unsigned en8[] = { 0,4, 12, 35,39 };
  static const unsigned en9[] = { 1,5, 13, 36,40 };
  static const unsigned en10[]= { 2,6, 14, 37,41 };
  static const unsigned en11[]= { 3,7, 15, 38,42 };

  static const unsigned * const edge_map[] =
    { en0, en1, en2, en3, en4, en5, en6, en7, en8, en9, en10, en11 };

  return edge_map ;
}

const unsigned * const * hex125_face_map()
{
  StaticAssert< 0 == Hex::side<0>::vertex<0>::ordinal >::ok();
  StaticAssert< 1 == Hex::side<0>::vertex<1>::ordinal >::ok();
  StaticAssert< 5 == Hex::side<0>::vertex<2>::ordinal >::ok();
  StaticAssert< 4 == Hex::side<0>::vertex<3>::ordinal >::ok();

  StaticAssert< 1 == Hex::side<1>::vertex<0>::ordinal >::ok();
  StaticAssert< 2 == Hex::side<1>::vertex<1>::ordinal >::ok();
  StaticAssert< 6 == Hex::side<1>::vertex<2>::ordinal >::ok();
  StaticAssert< 5 == Hex::side<1>::vertex<3>::ordinal >::ok();

  StaticAssert< 2 == Hex::side<2>::vertex<0>::ordinal >::ok();
  StaticAssert< 3 == Hex::side<2>::vertex<1>::ordinal >::ok();
  StaticAssert< 7 == Hex::side<2>::vertex<2>::ordinal >::ok();
  StaticAssert< 6 == Hex::side<2>::vertex<3>::ordinal >::ok();

  StaticAssert< 0 == Hex::side<3>::vertex<0>::ordinal >::ok();
  StaticAssert< 4 == Hex::side<3>::vertex<1>::ordinal >::ok();
  StaticAssert< 7 == Hex::side<3>::vertex<2>::ordinal >::ok();
  StaticAssert< 3 == Hex::side<3>::vertex<3>::ordinal >::ok();

  StaticAssert< 0 == Hex::side<4>::vertex<0>::ordinal >::ok();
  StaticAssert< 3 == Hex::side<4>::vertex<1>::ordinal >::ok();
  StaticAssert< 2 == Hex::side<4>::vertex<2>::ordinal >::ok();
  StaticAssert< 1 == Hex::side<4>::vertex<3>::ordinal >::ok();

  StaticAssert< 4 == Hex::side<5>::vertex<0>::ordinal >::ok();
  StaticAssert< 5 == Hex::side<5>::vertex<1>::ordinal >::ok();
  StaticAssert< 6 == Hex::side<5>::vertex<2>::ordinal >::ok();
  StaticAssert< 7 == Hex::side<5>::vertex<3>::ordinal >::ok();

  // Each face will have at 25 nodes to cover all cases of refined
  // Hexahedron faces: Quad_4_3D, Quad_8_3D, or Quad_9_3D.
  // Face node tables use Six groups of nodes:
  // a) vertices                 (Local nodes: 0-1-2-3               )
  // b) edge mid-points          (Local nodes: 4-5-6-7               )
  // c) centroid                 (Local node : 8                     )
  // d) edge quater points       (Local nodes: 9-10-11-12-13-14-15-16)
  // e) interior edge mid-points (Local nodes: 17-18-19-20           )
  // f) mid-quadrant points      (Local nodes: 21-22-23-24           )

  static const unsigned face_0[] = { 0,1,5,4,   8,13,16,12,  25,
                                 27,28,36,40,44,43,39,35,
                                 59,52,66,51,  105,106,107,108 };

  static const unsigned face_1[] = { 1,2,6,5,   9,14,17,13,  24,
                                 29,30,37,41,46,45,40,36,
                                 69,54,70,53,  101,102,103,104 };

  static const unsigned face_2[] = { 2,3,7,6,  10,15,18,14,  26,
                                 31,32,38,42,48,47,41,37,
                                 62,56,63,55,  109,110,111,112 };

  static const unsigned face_3[] = { 0,4,7,3, 12,19,15,11,  23,
                                 35,39,50,49,42,38,33,34,
                                 58,73,57,74, 97,98,99,100   };

  static const unsigned face_4[] = { 0,3,2,1,   11,10,9,8,   21,
                                 34,33,32,31,30,29,28,27,
                                 67,61,68,60,  89,90,91,92    };

  static const unsigned face_5[] = { 4,5,6,7,  16,17,18,19,  22,
                                 43,44,45,46,47,48,49,50,
                                 65,71,64,72,   93,94,95,96    };

  static const unsigned * const face_map[] =
    { face_0 , face_1 , face_2 , face_3 , face_4 , face_5 };

  return face_map ;
}


const Topology::Boundary * hex125_edge_2_map()
{
  static const Topology & edge_top = topology< Edge , 2 >();

  static const unsigned * const * const edge_map = hex125_edge_map();

  static Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top },
      { edge_map[3] , & edge_top },
      { edge_map[4] , & edge_top },
      { edge_map[5] , & edge_top },
      { edge_map[6] , & edge_top },
      { edge_map[7] , & edge_top },
      { edge_map[8] , & edge_top },
      { edge_map[9] , & edge_top },
      { edge_map[10] , & edge_top },
      { edge_map[11] , & edge_top } };

  return edges ;
}

const Topology::Boundary * hex125_edge_3_map()
{
  static const Topology & edge_top = topology< Edge , 3 >();

  static const unsigned * const * const edge_map = hex125_edge_map();

  static Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top },
      { edge_map[3] , & edge_top },
      { edge_map[4] , & edge_top },
      { edge_map[5] , & edge_top },
      { edge_map[6] , & edge_top },
      { edge_map[7] , & edge_top },
      { edge_map[8] , & edge_top },
      { edge_map[9] , & edge_top },
      { edge_map[10] , & edge_top },
      { edge_map[11] , & edge_top } };

  return edges ;
}

const Topology::Boundary * hex125_face_4_map()
{
  static const Topology & face_top = topology< QuadrilateralFace , 4 >();

  static const unsigned * const * const face_map = hex125_face_map();

  static Topology::Boundary faces[] =
    { { face_map[0] , & face_top } ,
      { face_map[1] , & face_top } ,
      { face_map[2] , & face_top } ,
      { face_map[3] , & face_top } ,
      { face_map[4] , & face_top } ,
      { face_map[5] , & face_top } };

  return faces ;
}

const Topology::Boundary * hex125_face_8_map()
{
  static const Topology & face_top = topology< QuadrilateralFace , 8 >();
  static const unsigned * const * const face_map = hex125_face_map();

  static Topology::Boundary faces[] =
    { { face_map[0] , & face_top } ,
      { face_map[1] , & face_top } ,
      { face_map[2] , & face_top } ,
      { face_map[3] , & face_top } ,
      { face_map[4] , & face_top } ,
      { face_map[5] , & face_top } };

  return faces ;
}

const Topology::Boundary * hex125_face_9_map()
{
  static const Topology & face_top = topology< QuadrilateralFace , 9 >();
  static const unsigned * const * const face_map = hex125_face_map();

  static Topology::Boundary faces[] =
    { { face_map[0] , & face_top } ,
      { face_map[1] , & face_top } ,
      { face_map[2] , & face_top } ,
      { face_map[3] , & face_top } ,
      { face_map[4] , & face_top } ,
      { face_map[5] , & face_top } };

  return faces ;
}

}

//----------------------------------------------------------------------

template<>
const Topology & topology< Hex , 8 >()
{
  static const Topology top( Hex(), 8 , 125 , hex125_edge_2_map(),
                                              hex125_face_4_map() );
  return top ;
}

template<>
const Topology & topology< Hex , 20 >()
{
  static const Topology top( Hex(), 20 , 125 , hex125_edge_3_map(),
                                               hex125_face_8_map() );

  return top ;
}

template<>
const Topology & topology< Hex , 27 >()
{
  static const Topology top( Hex(), 20 , 125 , hex125_edge_3_map(),
                                               hex125_face_9_map() );
  return top ;
}

//----------------------------------------------------------------------

}
}

