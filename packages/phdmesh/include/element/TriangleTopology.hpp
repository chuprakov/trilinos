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

#ifndef phdmesh_TriangleTopology_hpp
#define phdmesh_TriangleTopology_hpp

#include <element/LocalTopology.hpp>

namespace phdmesh {

template< unsigned NumNode = 0 > struct BoundaryTriangle ; // 3D only
template< unsigned NumNode = 0 > struct Triangle ;         // 2D only
template< unsigned NumNode = 0 > struct TriangleShell ;    // 3D only

/*----------------------------------------------------------------------*/
/*  Triangle node numbering, up to 6 nodes:                             */
/*                                                                      */
/*           2                                                          */
/*           o                                                          */
/*          / \                                                         */
/*         /   \                                                        */
/*        /     \                                                       */
/*     5 o       o 4                                                    */
/*      /         \                                                     */
/*     /           \                                                    */
/*    /             \                                                   */
/*   o-------o-------o                                                  */
/*  0        3        1                                                 */
/*                                                                      */
/*  Triangle node numbering, up to 15 nodes conformal to 6 node:        */
/*                                                                      */
/*           2                                                          */
/*           o                                                          */
/*          / \                                                         */
/*      10 *   * 9                                                      */
/*        / 14  \                                                       */
/*     5 o---*---o 4                                                    */
/*      / \     / \                                                     */
/*  11 * 12*   *13 * 8                                                  */
/*    /     \ /     \                                                   */
/*   o---*---o---*---o                                                  */
/*  0    6   3   7    1                                                 */
/*                                                                      */
/*  Interior node 12 is on the line between node 0 and node 4           */
/*  Interior node 13 is on the line between node 1 and node 5           */
/*  Interior node 14 is on the line between node 2 and node 3           */
/*                                                                      */
/*----------------------------------------------------------------------*/

typedef
  MakeTypeList< IndexList< 0 , 1 , 3 ,  6 ,  7 > ,
                IndexList< 1 , 2 , 4 ,  8 ,  9 > ,
                IndexList< 2 , 0 , 5 , 10 , 11 > >::type TriangleEdgeNodeMap ;

template<>
struct BoundaryTriangle<0>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                3 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      TriangleEdgeNodeMap >
{
  typedef BoundaryTriangle<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryTriangle<3>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                3 , 3, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      TriangleEdgeNodeMap >
{
  typedef BoundaryTriangle<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryTriangle<6>
  : public LocalTopologyTraits< 2 , 3 , 3, // Dimension
                                3 , 6, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      TriangleEdgeNodeMap >
{
  typedef BoundaryTriangle<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------

template<>
struct Triangle<0>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                3 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      TriangleEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      TriangleEdgeNodeMap >
{
  typedef Triangle<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Triangle<3>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                3 , 3, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      TriangleEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      TriangleEdgeNodeMap >
{
  typedef Triangle<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Triangle<6>
  : public LocalTopologyTraits< 2 , 2 , 2, // Dimension
                                3 , 6, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      TriangleEdgeNodeMap ,
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      TriangleEdgeNodeMap >
{
  typedef Triangle<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 , 2 > ,
                IndexList< 0 , 2 , 1 > >::type TriangleShellSideNodeMap ;

template<>
struct TriangleShell<0>
  : public LocalTopologyTraits< 3 , 3 , 3, // Dimension
                                3 , 0, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      TriangleEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<> , BoundaryTriangle<> >::type ,
      TriangleShellSideNodeMap >
{
  typedef TriangleShell<0> Shape ; 
  static const LocalTopology * descriptor();
};

template<>
struct TriangleShell<3>
  : public LocalTopologyTraits< 3 , 3 , 3, // Dimension
                                3 , 3, // #Vertex, #Node
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      TriangleEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<3> , BoundaryTriangle<3> >::type ,
      TriangleShellSideNodeMap >
{
  typedef TriangleShell<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

typedef MakeTypeList<
  /* Rotation = 0 , Polarity = true */
  IndexList< 0 , 1 , 2 , 3 , 4 , 5 ,
             6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 > ,
  /* Rotation = 1 , Polarity = true */
  IndexList< 2 ,  0 ,  1 ,  5 ,  3 ,  4 ,
            10 , 11 ,  6 ,  7 ,  8 ,  9 , 14 , 12 , 13 > ,
  /* Rotation = 2 , Polarity = true */
  IndexList< 1 ,  2 ,  0 , 4 ,  5 ,  3 ,
             8 ,  9 , 10 , 11 ,  6 ,  7 , 13 , 14 , 12 > ,
  /* Rotation = 0 , Polarity = false */
  IndexList< 0 ,  2 , 1 , 5 , 4 , 3 ,
            11 , 10 , 9 , 8 , 7 , 6 , 12 , 14 , 13 > ,
  /* Rotation = 1 , Polarity = false */
  IndexList< 2 ,  1 ,  0 ,  4 ,  3 ,  5 ,
             9 ,  8 ,  7 ,  6 , 11 , 10 , 14 , 13 , 12 > ,
  /* Rotation = 2 , Polarity = false */
  IndexList< 1 ,  0 ,  2 , 3 ,  5 ,  4 ,
             7 ,  6 , 11 , 10 ,  9 ,  8 , 13 , 12 , 14 >
  >::type TriangleRotationNodeMap ;

} // namespace phdmesh


#endif

