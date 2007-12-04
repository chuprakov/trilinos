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

typedef TriangleFace TFace ;

namespace {

const unsigned * const * tri15_edge_map()
{
  StaticAssert< 3 == TFace::number_vertex >::ok();
  StaticAssert< 3 == TFace::number_edge >::ok();
  StaticAssert< 0 == TFace::number_side >::ok();
  StaticAssert<      TFace::boundary >::ok();
  StaticAssert< 0 == TFace::edge<0>::vertex<0>::ordinal >::ok();
  StaticAssert< 1 == TFace::edge<0>::vertex<1>::ordinal >::ok();
  StaticAssert< 1 == TFace::edge<1>::vertex<0>::ordinal >::ok();
  StaticAssert< 2 == TFace::edge<1>::vertex<1>::ordinal >::ok();
  StaticAssert< 2 == TFace::edge<2>::vertex<0>::ordinal >::ok();
  StaticAssert< 0 == TFace::edge<2>::vertex<1>::ordinal >::ok();

  static const unsigned en0[] = { 0, 1,  3,   6,  7 };
  static const unsigned en1[] = { 1, 2,  4,   8,  9 };
  static const unsigned en2[] = { 2, 0,  5,  10, 11 };
  static const unsigned * const tri_map[] = { en0 , en1 , en2 };
  return tri_map ;
}

const Topology::Boundary * tri15_edge_2_map()
{
  static const Topology & edge_top = topology< Edge , 2 >();

  static const unsigned * const * const edge_map = tri15_edge_map();

  static const Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top } };

  return edges ;
}

const Topology::Boundary * tri15_edge_3_map()
{
  static const Topology & edge_top = topology< Edge , 3 >();

  static const unsigned * const * const edge_map = tri15_edge_map();

  static const Topology::Boundary edges[] =
    { { edge_map[0] , & edge_top },
      { edge_map[1] , & edge_top },
      { edge_map[2] , & edge_top } };

  return edges ;
}

const Topology::Orientation * tri15_orientation()
{
  StaticAssert< 6 == Orientation< TFace >::number_orientation >::ok();

  StaticAssert< 0 == Orientation< TFace , 0 , true  >::ordinal >::ok() ;
  StaticAssert< 1 == Orientation< TFace , 1 , true  >::ordinal >::ok() ;
  StaticAssert< 2 == Orientation< TFace , 2 , true  >::ordinal >::ok() ;
  StaticAssert< 3 == Orientation< TFace , 0 , false >::ordinal >::ok() ;
  StaticAssert< 4 == Orientation< TFace , 1 , false >::ordinal >::ok() ;
  StaticAssert< 5 == Orientation< TFace , 2 , false >::ordinal >::ok() ;

  static const unsigned tri15_R0_DP[] = { 0 , 1 , 2 ,
                                          3 , 4 , 5 ,
                                          6 , 7 , 8 , 9 , 10 , 11 ,
                                         12 , 13 , 14 };

  static const unsigned tri15_R1_DP[] = { 2 ,  0 ,  1 ,
                                          5 ,  3 ,  4 ,
                                         10 , 11 ,  6 ,  7 ,  8 ,  9 ,
                                         14 , 12 , 13 };

  static const unsigned tri15_R2_DP[] = { 1 ,  2 ,  0 ,
                                          4 ,  5 ,  3 ,
                                          8 ,  9 , 10 , 11 ,  6 ,  7 ,
                                         13 , 14 , 12 };

  static const unsigned tri15_R0_DM[] = { 0 , 2 , 1 ,
                                          5 , 4 , 3 ,
                                         11 , 10 , 9 , 8 , 7 , 6 ,
                                         12 , 14 , 13 };

  static const unsigned tri15_R1_DM[] = { 2 ,  1 ,  0 ,
                                          4 ,  3 ,  5 ,
                                          9 ,  8 ,  7 ,  6 , 11 , 10 ,
                                         14 , 13 , 12 };

  static const unsigned tri15_R2_DM[] = { 1 ,  0 ,  2 ,
                                          3 ,  5 ,  4 ,
                                          7 ,  6 , 11 , 10 ,  9 ,  8 ,
                                         13 , 12 , 14 };

  static Topology::Orientation orientation[] =
    { { tri15_R0_DP , 0 , true },
      { tri15_R1_DP , 1 , true },
      { tri15_R2_DP , 2 , true },
      { tri15_R0_DM , 0 , false },
      { tri15_R1_DM , 1 , false },
      { tri15_R2_DM , 2 , false } };

  return orientation ;
}

}

//----------------------------------------------------------------------

template<>
const Topology & topology< TFace , 3 >()
{
  static const Topology top( TFace() , 3 , 15 ,
                             tri15_edge_2_map(), NULL ,
                             tri15_orientation() );
  return top ;
}

template<>
const Topology & topology< TFace , 6 >()
{
  static const Topology top( TFace() , 6 , 15 ,
                             tri15_edge_3_map(), NULL ,
                             tri15_orientation() );
  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< Triangle2D , 3 >()
{
  static const Topology top( Triangle2D() , 3 , 15 , tri15_edge_2_map() );
  return top ;
}

template<>
const Topology & topology< Triangle2D , 6 >()
{
  static const Topology top( Triangle2D() , 6 , 15 , tri15_edge_3_map() );
  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< TriangleShell , 3 >()
{
  static const Topology & face_top = topology< TFace , 3 >();

  static Topology::Boundary faces[] =
    { { face_top.orientation[0].node_map , & face_top } ,
      { face_top.orientation[3].node_map , & face_top } };

  static const Topology top( TriangleShell(), 3, 15,
                             tri15_edge_2_map() , faces );

  return top ;
}

//----------------------------------------------------------------------

template<>
const Topology & topology< TriangleShell , 6 >()
{
  static const Topology & face_top = topology< TFace , 6 >();

  static Topology::Boundary faces[] =
    { { face_top.orientation[0].node_map , & face_top } ,
      { face_top.orientation[3].node_map , & face_top } };

  static const Topology top( TriangleShell(), 6, 15,
                             tri15_edge_3_map() , faces );

  return top ;
}

//----------------------------------------------------------------------

}
}

