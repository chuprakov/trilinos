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

#ifndef phdmesh_PyramidTopology_hpp
#define phdmesh_PyramidTopology_hpp

#include <element/LocalTopology.hpp>
#include <element/TriangleTopology.hpp>
#include <element/QuadrilateralTopology.hpp>

namespace phdmesh {

template< unsigned NumNode = 0 > struct Pyramid ;      // 3D only

//----------------------------------------------------------------------

typedef
  MakeTypeList< IndexList< 0 , 1 ,   5 > ,
                IndexList< 1 , 2 ,   6 > ,
                IndexList< 2 , 0 ,   7 > ,
                IndexList< 0 , 3 ,   8 > ,
                IndexList< 0 , 4 ,   9 > ,
                IndexList< 1 , 4 ,  10 > ,
                IndexList< 2 , 4 ,  11 > ,
                IndexList< 3 , 4 ,  12 > >::type PyramidEdgeNodeMap ;

typedef
  MakeTypeList< IndexList< 0 , 1 , 4 ,       5 , 10 ,  9 > ,
                IndexList< 1 , 2 , 4 ,       6 , 11 , 10 > ,
                IndexList< 2 , 3 , 4 ,       7 , 12 , 11 > ,
                IndexList< 3 , 0 , 4 ,       8 ,  9 , 12 > ,
                IndexList< 0 , 3 , 2 , 1 ,   8 ,  7 ,  6 ,  5 ,  13 >
  >::type PyramidFaceNodeMap ;

template<>
struct Pyramid<0>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                5 , 0 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> ,
                    BoundaryEdge<> >::type ,
      PyramidEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<> ,
                    BoundaryTriangle<> ,
                    BoundaryTriangle<> ,
                    BoundaryTriangle<> ,
                    BoundaryQuadrilateral<> >::type ,
      PyramidFaceNodeMap >
{
  typedef Pyramid<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Pyramid<5>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                5 , 5 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> ,
                    BoundaryEdge<2> >::type ,
      PyramidEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<3> ,
                    BoundaryTriangle<3> ,
                    BoundaryTriangle<3> ,
                    BoundaryTriangle<3> ,
                    BoundaryQuadrilateral<4> >::type ,
      PyramidFaceNodeMap >
{
  typedef Pyramid<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Pyramid<13>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                5 , 13 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      PyramidEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryQuadrilateral<8> >::type ,
      PyramidFaceNodeMap >
{
  typedef Pyramid<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Pyramid<14>
  : public LocalTopologyTraits< 3 , 3 , 3 , // Dimension
                                5 , 14 , // #Vertex, #Node ,
      MakeTypeList< BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> ,
                    BoundaryEdge<3> >::type ,
      PyramidEdgeNodeMap ,
      MakeTypeList< BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryTriangle<6> ,
                    BoundaryQuadrilateral<9> >::type ,
      PyramidFaceNodeMap >
{
  typedef Pyramid<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

