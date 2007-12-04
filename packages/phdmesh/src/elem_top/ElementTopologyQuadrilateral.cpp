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
// Conventional node numbering for nine-node quadrilateral
//
//    3        6        2
//     o-------o-------o
//     |               |
//     |               |
//     |       8       |
//   7 o       o       o 5
//     |               |
//     |               |
//     |               |
//     o-------o-------o
//    0        4        1
//
//
// Conformal node numbering with up to 25 nodes:
//
//    3    14    6   13     2
//     o----*----o----*----o
//     |         |         |
//     |   24    |    23   |
//   15*    *    *19  *    *12
//     |         |         |
//     |        8|    18   |
//   7 o----*----o----*----o 5
//     |   20    |         |
//     |         |         |
//   16*    *  17*    *    *11
//     |   21    |   22    |
//     |         |         |
//     o----*----o----*----o
//    0     9    4   10     1
//
//----------------------------------------------------------------------
// Conventional edge numbering for quadrilateral
//
//                 Edge #2
//           3        6        2
//            o-------o-------o
//            |               |
//            |               |
//            |       8       |
// Edge #3  7 o       o       o 5  Edge #1
//            |               |
//            |               |
//            |               |
//            o-------o-------o
//           0        4        1
//                  Edge #0
//
//----------------------------------------------------------------------

typedef QuadrilateralFace QFace ;

namespace {

const unsigned * const * quad25_edge_map()
{
  StaticAssert< 4 == QFace::number_vertex >::ok();
  StaticAssert< 4 == QFace::number_edge >::ok();
  StaticAssert< 0 == QFace::number_side >::ok();
  StaticAssert<      QFace::boundary >::ok();
  StaticAssert< 0 == QFace::edge<0>::vertex<0>::ordinal >::ok();
  StaticAssert< 1 == QFace::edge<0>::vertex<1>::ordinal >::ok();
  StaticAssert< 1 == QFace::edge<1>::vertex<0>::ordinal >::ok();
  StaticAssert< 2 == QFace::edge<1>::vertex<1>::ordinal >::ok();
  StaticAssert< 2 == QFace::edge<2>::vertex<0>::ordinal >::ok();
  StaticAssert< 3 == QFace::edge<2>::vertex<1>::ordinal >::ok();
  StaticAssert< 3 == QFace::edge<3>::vertex<0>::ordinal >::ok();
  StaticAssert< 0 == QFace::edge<3>::vertex<1>::ordinal >::ok();

  static const unsigned en0[] = { 0, 1,  4,   9, 10 };
  static const unsigned en1[] = { 1, 2,  5,  11, 12 };
  static const unsigned en2[] = { 2, 3,  6,  13, 14 };
  static const unsigned en3[] = { 3, 0,  7,  15, 16 };

  static const unsigned * const edge_map[] = { en0 , en1 , en2 , en3 };

  return edge_map ;
}

const Topology::Boundary * quad25_edge_2_map()
{
  static const Topology & edge_top = topology< Edge , 2 >();

  static const unsigned * const * const edge_map = quad25_edge_map();

  static Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top },
      { edge_map[3] , & edge_top } };

  return edges ;
}

const Topology::Boundary * quad25_edge_3_map()
{
  static const Topology & edge_top = topology< Edge , 3 >();

  static const unsigned * const * const edge_map = quad25_edge_map();

  static Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top },
      { edge_map[3] , & edge_top } };

  return edges ;
}

const Topology::Orientation * quad25_orientation()
{
  StaticAssert< 8 == Orientation<QFace>::number_orientation >::ok();

  StaticAssert< 0 == Orientation<QFace, 0 , true  >::ordinal >::ok() ;
  StaticAssert< 1 == Orientation<QFace, 1 , true  >::ordinal >::ok() ;
  StaticAssert< 2 == Orientation<QFace, 2 , true  >::ordinal >::ok() ;
  StaticAssert< 3 == Orientation<QFace, 3 , true  >::ordinal >::ok() ;
  StaticAssert< 4 == Orientation<QFace, 0 , false >::ordinal >::ok() ;
  StaticAssert< 5 == Orientation<QFace, 1 , false >::ordinal >::ok() ;
  StaticAssert< 6 == Orientation<QFace, 2 , false >::ordinal >::ok() ;
  StaticAssert< 7 == Orientation<QFace, 3 , false >::ordinal >::ok() ;

  static const unsigned np_R0_DP[] = {  0,    1,    2,    3,
                                        4,    5,    6,    7,    8,
                                       9,10,11,12,13,14,15,16,
                                       17,   18,   19,   20,
                                       21,   22,   23,   24 };

  static const unsigned np_R0_DM[] = {  0,    3,    2,    1,
                                        7,    6,    5,    4,    8,
                                      16,15,14,13,12,11,10, 9,
                                       20,   19,   18,   17,
                                       21,   24,   23,   22 };

  static const unsigned np_R1_DP[] = { 3,    0,    1,    2,
                                       7,    4,    5,    6,    8,
                                     15,16, 9,10,11,12,13,14,
                                      20,   17,   18,   19,
                                      24,   21,   22,   23 };

  static const unsigned np_R1_DM[] = { 3,    2,    1,    0,
                                       6,    5,    4,    7,    8,
                                     14,13,12,11,10, 9,16,15,
                                      19,   18,   17,   20,
                                      24,   23,   22,   21 };

  static const unsigned np_R2_DP[] = { 2,    3,    0,    1,
                                       6,    7,    4,    5,    8,
                                     13,14,15,16, 9,10,11,12,
                                      19,   20,   17,   18,
                                      23,   24,   21,   22 };

  static const unsigned np_R2_DM[] = { 2,    1,    0,    3,
                                       5,    4,    7,    6,    8,
                                     12,11,10, 9,16,15,14,13,
                                      18,   17,   20,   19,
                                      23,   22,   21,   24 };

  static const unsigned np_R3_DP[] = { 1,    2,    3,    0,
                                       5,    6,    7,    4,    8,
                                     11,12,13,14,15,16, 9,10,
                                      18,   19,   20,   17,
                                      22,   23,   24,   21 };

  static const unsigned np_R3_DM[] = { 1,    0,    3,    2,
                                       4,    7,    6,    5,    8,
                                     10, 9,16,15,14,13,12,11,
                                      17,   20,   19,   18,
                                      22,   21,   24,   23 };

  static const Topology::Orientation orientation[] =
    { { np_R0_DP , 0 , true },
      { np_R1_DP , 1 , true },
      { np_R2_DP , 2 , true },
      { np_R3_DP , 3 , true },
      { np_R0_DM , 0 , false },
      { np_R1_DM , 1 , false },
      { np_R2_DM , 2 , false },
      { np_R3_DM , 3 , false } };

  return orientation ;
}

}

//----------------------------------------------------------------------

template<>
const Topology & topology< QFace , 4 >()
{
  static const Topology
    top( QFace(), 4 , 25 , quad25_edge_2_map(), NULL, quad25_orientation() );

  return top ;
}

template<>
const Topology & topology< QFace , 8 >()
{
  static const Topology
    top( QFace(), 8 , 25 , quad25_edge_3_map(), NULL, quad25_orientation() );

  return top ;
}

template<>
const Topology & topology< QFace , 9 >()
{
  static const Topology
    top( QFace(), 9 , 25 , quad25_edge_3_map(), NULL, quad25_orientation() );

  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< Quadrilateral2D , 3 >()
{
  static const Topology top( Quadrilateral2D(), 4, 25, quad25_edge_2_map() );
  return top ;
}

template<>
const Topology & topology< Quadrilateral2D , 8 >()
{
  static const Topology top( Quadrilateral2D(), 8, 25, quad25_edge_3_map() );
  return top ;
}

template<>
const Topology & topology< Quadrilateral2D , 9 >()
{
  static const Topology top( Quadrilateral2D(), 9, 25, quad25_edge_3_map() );
  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< QuadrilateralShell , 4 >()
{
  static const Topology & face_top = topology< QFace , 4 >();

  static Topology::Boundary faces[] =
    { { face_top.orientation[0].node_map , & face_top } ,
      { face_top.orientation[3].node_map , & face_top } };

  static const Topology top( QuadrilateralShell(), 4, 25,
                             quad25_edge_2_map(), faces );

  return top ;
}

template<>
const Topology & topology< QuadrilateralShell , 8 >()
{
  static const Topology & face_top = topology< QFace , 8 >();

  static Topology::Boundary faces[] =
    { { face_top.orientation[0].node_map , & face_top } ,
      { face_top.orientation[3].node_map , & face_top } };

  static const Topology top( QuadrilateralShell(), 8, 25,
                             quad25_edge_3_map(), faces );

  return top ;
}

template<>
const Topology & topology< QuadrilateralShell , 9 >()
{
  static const Topology & face_top = topology< QFace , 9 >();

  static Topology::Boundary faces[] =
    { { face_top.orientation[0].node_map , & face_top } ,
      { face_top.orientation[3].node_map , & face_top } };

  static const Topology top( QuadrilateralShell(), 9, 25,
                             quad25_edge_3_map(), faces );

  return top ;
}

//----------------------------------------------------------------------

}
}

