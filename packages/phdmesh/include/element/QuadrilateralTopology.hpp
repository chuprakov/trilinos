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

#ifndef phdmesh_QuadrilateralTopology_hpp
#define phdmesh_QuadrilateralTopology_hpp

#include <element/LocalTopology.hpp>

namespace phdmesh {

template< unsigned NumNode = 0 > struct BoundaryQuadrilateral ; // 3D only
template< unsigned NumNode = 0 > struct Quadrilateral ;         // 2D only
template< unsigned NumNode = 0 > struct QuadrilateralShell ;    // 3D only

//----------------------------------------------------------------------
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
//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 ,  9 , 10 > ,
                IndexList< 1 , 2 , 5 , 11 , 12 > ,
                IndexList< 2 , 3 , 6 , 13 , 14 > ,
                IndexList< 3 , 0 , 7 , 15 , 16 > >::type
  QuadrilateralEdgeNodeMap ;

template<>
struct BoundaryQuadrilateral<0>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                4 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef BoundaryQuadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryQuadrilateral<4>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                4 , 4, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef BoundaryQuadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryQuadrilateral<8>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                4 , 8, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef BoundaryQuadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryQuadrilateral<9>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                4 , 9, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef BoundaryQuadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------

template<>
struct Quadrilateral<0>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                4 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> , BoundaryEdge<> ,
                    BoundaryEdge<> , BoundaryEdge<> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<> , BoundaryEdge<> ,
                    BoundaryEdge<> , BoundaryEdge<> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef Quadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Quadrilateral<4>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                4 , 4, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> , BoundaryEdge<2> ,
                    BoundaryEdge<2> , BoundaryEdge<2> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<2> , BoundaryEdge<2> ,
                    BoundaryEdge<2> , BoundaryEdge<2> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef Quadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Quadrilateral<8>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                4 , 8, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> , BoundaryEdge<3> ,
                    BoundaryEdge<3> , BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<3> , BoundaryEdge<3> ,
                    BoundaryEdge<3> , BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef Quadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Quadrilateral<9>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                4 , 9, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> , BoundaryEdge<3> ,
                    BoundaryEdge<3> , BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<3> , BoundaryEdge<3> ,
                    BoundaryEdge<3> , BoundaryEdge<3> >::type ,
      QuadrilateralEdgeNodeMap >
{
  typedef Quadrilateral<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 2 , 3 > ,
                IndexList< 0 , 3 , 2 , 1 > >::type
    QuadrilateralShellSideNodeMap ;

template<>
struct QuadrilateralShell<0>
  : public LocalTopologyTraits< 3 , 3 , 3, // Dimension
                                4 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<> ,
                    BoundaryQuadrilateral<> >::type ,
      QuadrilateralShellSideNodeMap >
{
  typedef QuadrilateralShell<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct QuadrilateralShell<4>
  : public LocalTopologyTraits< 3 , 3 , 3, // Dimension
                                4 , 4, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      QuadrilateralEdgeNodeMap ,
      MakeTypeList< BoundaryQuadrilateral<4> ,
                    BoundaryQuadrilateral<4> >::type ,
      QuadrilateralShellSideNodeMap >
{
  typedef QuadrilateralShell<0> Shape ; 
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

typedef
  MakeTypeList<
  /* Rotation = 0 , Polarity = true */
  IndexList<  0,  1,  2,  3,  4,  5,  6,  7,  8,
              9, 10, 11, 12, 13, 14, 15, 16,
             17, 18, 19, 20, 21, 22, 23, 24 > ,
  /* Rotation = 1 , Polarity = true */
  IndexList<  3,  0,  1,  2,  7,  4,  5,  6,  8,
             15, 16,  9, 10, 11, 12, 13, 14,
             20, 17, 18, 19, 24, 21, 22, 23 > ,
  /* Rotation = 2 , Polarity = true */
  IndexList<  2,  3,  0,  1,  6,  7,  4,  5,  8,
             13, 14, 15, 16,  9, 10, 11, 12,
             19, 20, 17, 18, 23, 24, 21, 22 > ,
  /* Rotation = 3 , Polarity = true */
  IndexList<  1,  2,  3,  0,  5,  6,  7,  4,  8,
             11, 12, 13, 14, 15, 16,  9, 10,
             18, 19, 20, 17, 22, 23, 24, 21 > ,
  /* Rotation = 0 , Polarity = false */
  IndexList<  0,  3,  2,  1,  7,  6,  5,  4,  8,
             16, 15, 14, 13, 12, 11, 10,  9,
             20, 19, 18, 17, 21, 24, 23, 22 > ,
  /* Rotation = 1 , Polarity = false */
  IndexList<  3,  2,  1,  0,  6,  5,  4,  7,  8,
             14, 13, 12, 11, 10,  9, 16, 15,
             19, 18, 17, 20, 24, 23, 22, 21 > ,
  /* Rotation = 2 , Polarity = false */
  IndexList<  2,  1,  0,  3,  5,  4,  7,  6,  8,
             12, 11, 10,  9, 16, 15, 14, 13,
             18, 17, 20, 19, 23, 22, 21, 24 > ,
  /* Rotation = 3 , Polarity = false */
  IndexList<  1,  0,  3,  2,  4,  7,  6,  5,  8,
             10,  9, 16, 15, 14, 13, 12, 11,
             17, 20, 19, 18, 22, 21, 24, 23 >
  >::type QuadrilateralRotationNodeMap ;

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

