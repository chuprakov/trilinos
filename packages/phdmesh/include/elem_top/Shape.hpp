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

#ifndef phdmesh_ElementShape_hpp
#define phdmesh_ElementShape_hpp

#include <util/TypeList.hpp>

namespace phdmesh {
namespace element {

//----------------------------------------------------------------------
/** @class Shape
 *  @brief Traits for an element shape.
 *
 *  Shape::number_vertex
 *  Shape::number_edge
 *  Shape::number_side                  duplicates edges in 2D
 *  Shape::boundary                     if describing an element boundary
 *  Shape::mimimum_dimension            minimum valid spatial dimension
 *  Shape::maximum_dimension            maximum valid spatial dimension
 * 
 *  Shape::edge<I>                      defined for I < number_edge
 *  Shape::edge<I>::number_vertex
 *  Shape::edge<I>::number_edge
 *  Shape::edge<I>::number_side
 *  Shape::edge<I>::boundary            == true
 *  Shape::edge<I>::vertex<J>::ordinal  defined for J < edge<I>::number_vertex
 * 
 *  Shape::side<I>                      defined for I < number_side
 *  Shape::side<I>::number_vertex
 *  Shape::side<I>::number_edge
 *  Shape::side<I>::number_side
 *  Shape::side<I>::boundary            == true
 *  Shape::side<I>::vertex<J>::ordinal  defined for J < side<I>::number_vertex
 */
template< unsigned Number_Vertex ,
          class EdgeList ,
          class SideList ,
          unsigned MinimumDimension ,
          unsigned MaximumDimension ,
          bool Boundary = false > struct Shape ;

//----------------------------------------------------------------------
/** @class Orientation
 *  @brief Potential orientations of an element boundary shape.
 *
 *  Orientation< Shape >::number_orientation
 *  Orientation< BoundaryShape , Ordinal >::rotation
 *  Orientation< BoundaryShape , Ordinal >::polarity
 *  Orientation< BoundaryShape , Rotation , Polarity >::ordinal
 *
 *  Polarity indicates whether the boundary has the expected direction,
 *  e.g. outward face normal.
 *
 *  Rotation indicates at which vertex the expected 0-th vertex
 *  actually appears.
 *
 *  Boundary shape orientation:
 *   nv = number of vertices
 *   jv = vertex appearing in the 'iv' location
 *   jv = polarity ? ( ( iv      ) + ( nv - rotation ) ) % nv
 *                 : ( ( nv - iv ) + ( nv - rotation ) ) % nv
 */
template< class S , unsigned Ordinal = 0 , bool Polarity = true >
class Orientation ;

//----------------------------------------------------------------------

template<class M, unsigned I> struct ShapeBoundaryVertex ;

template< unsigned Number_Vertex ,
          class EdgeList , // TypeList of ShapeBoundary types
          class SideList , // TypeList of ShapeBoundary types
          unsigned MinimumDimension ,
          unsigned MaximumDimension ,
          bool Boundary >
struct Shape {
  typedef Shape< Number_Vertex, EdgeList, SideList,
                 MinimumDimension, MaximumDimension,
                 Boundary> Self ;

  enum { number_vertex = Number_Vertex };
  enum { number_edge   = TypeListLength< EdgeList >::value };
  enum { number_side   = TypeListLength< SideList >::value };
  enum { boundary      = Boundary };

  enum { minimum_dimension = MinimumDimension };
  enum { maximum_dimension = MaximumDimension };

  template< unsigned I >
  class edge : public TypeListAt< EdgeList , I >::type {
    enum { OK = StaticAssert< I < Shape::number_edge >::OK };
  };

  template< unsigned I >
  class side : public TypeListAt< SideList , I >::type {
    enum { OK = StaticAssert< I < Shape::number_side >::OK };
  };

  ~Shape() {}
  Shape() {}
  Shape( const Self & ) {}
  Self & operator = ( const Self & ) {}
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Helper for defining boundary element relationships

template< unsigned J, unsigned v0, unsigned v1, unsigned v2, unsigned v3>
struct Select ;

template<unsigned v0, unsigned v1, unsigned v2, unsigned v3>
struct Select<0,v0,v1,v2,v3> { enum { value = v0 }; };

template<unsigned J, unsigned v0, unsigned v1, unsigned v2, unsigned v3>
struct Select { enum { value = Select<J-1,v1,v2,v3,0>::value }; };

template< class S , unsigned v0 = 0 , unsigned v1 = 0 ,
                    unsigned v2 = 0 , unsigned v3 = 0 >
struct ShapeBoundary : public S {
  template< unsigned J > struct vertex {
    enum { OK = StaticAssert< J < S::number_vertex >::OK };
    enum { ordinal = Select<J,v0,v1,v2,v3>::value };
  };
};

//----------------------------------------------------------------------
// Boundary shapes:

typedef Shape< 0 , TypeListEnd , TypeListEnd , 1 , 3 , true > Point ;
typedef Shape< 2 , TypeListEnd , TypeListEnd , 2 , 3 , true > Edge ;

typedef TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 0 > ,
        TypeListEnd > > >
        TriangleEdges ;

typedef Shape< 3 , TriangleEdges , TypeListEnd , 3 , 3 , true >
        TriangleFace ;

typedef TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 3 > ,
        TypeList< ShapeBoundary< Edge , 3 , 0 > ,
        TypeListEnd > > > >
        QuadrilateralEdges ;

typedef Shape< 4, QuadrilateralEdges, TypeListEnd, 3 , 3 , true >
        QuadrilateralFace ;

//----------------------------------------------------------------------

template< unsigned Ordinal , bool Polarity >
class Orientation<Edge,Ordinal,Polarity> {
private:
  enum { NV = Edge::number_vertex };
  enum { OK_Ordinal = StaticAssert< Ordinal < NV >::OK };
public:
  enum { number_orientation = 2 };
  enum { rotation = 1 == Ordinal || ! Polarity ? 1 : 0 };
  enum { polarity = 0 == rotation };
  enum { ordinal  = rotation };

  template<unsigned I>
  class vertex {
  private:
    enum { OK_I = StaticAssert< I < NV >::OK };
  public:
    enum { ordinal = polarity ? I : ( I + 1 ) % NV };
  };
};

template<unsigned Ordinal , bool Polarity >
class Orientation<TriangleFace,Ordinal,Polarity> {
private:
  enum { NV = TriangleFace::number_vertex };
  enum { OK_Ordinal = StaticAssert< Ordinal < NV * 2 >::OK };
public:
  enum { number_orientation = NV * 2 };
  enum { rotation = Ordinal % NV };
  enum { polarity = Ordinal < NV ? Polarity : false };
  enum { ordinal  = rotation + ( ! polarity ? NV : 0 ) };

  template<unsigned I>
  class vertex {
  private:
    enum { OK_I = StaticAssert< I < NV >::OK };
  public:
    enum { ordinal = ( ( NV - rotation ) + polarity ? I : NV - I ) % NV };
  };
};

template<unsigned Ordinal , bool Polarity >
class Orientation<QuadrilateralFace,Ordinal,Polarity> {
private:
  enum { NV = QuadrilateralFace::number_vertex };
  enum { OK_Ordinal = StaticAssert< Ordinal < NV * 2 >::OK };
public:
  enum { number_orientation = NV * 2 };
  enum { rotation = Ordinal % NV };
  enum { polarity = Ordinal < NV ? Polarity : false };
  enum { ordinal  = rotation + ( ! polarity ? NV : 0 ) };

  template<unsigned I>
  class vertex {
  private:
    enum { OK_I = StaticAssert< I < NV >::OK };
  public:
    enum { ordinal = ( ( NV - rotation ) + polarity ? I : NV - I ) % NV };
  };
};

template<class S, unsigned Ordinal, bool Polarity>
struct Orientation {
  enum { number_orientation = 0 };
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Solid element shapes:

typedef Shape< 1 , TypeListEnd , TypeListEnd , 1 , 3 >
        Sphere ;

//----------------------------------

typedef Shape< 2 , TypeListEnd , TypeListEnd , 1 , 1 >
        Line1D ;

//----------------------------------

typedef Shape< 3, TriangleEdges , TriangleEdges , 2 , 2 >
        Triangle2D ;

//----------------------------------

typedef Shape< 4, QuadrilateralEdges, QuadrilateralEdges , 2 , 2 >
        Quadrilateral2D ;  

//----------------------------------

typedef Shape< 4 ,
        // Edges:
        TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 0 > ,
        TypeList< ShapeBoundary< Edge , 0 , 3 > ,
        TypeList< ShapeBoundary< Edge , 1 , 3 > ,
        TypeList< ShapeBoundary< Edge , 2 , 3 > ,
        TypeListEnd > > > > > > ,
        // Faces:
        TypeList< ShapeBoundary< TriangleFace , 0 , 1 , 3 > ,
        TypeList< ShapeBoundary< TriangleFace , 1 , 2 , 3 > ,
        TypeList< ShapeBoundary< TriangleFace , 0 , 3 , 2 > ,
        TypeList< ShapeBoundary< TriangleFace , 0 , 2 , 1 > ,
        TypeListEnd > > > > ,
        3 , 3 > Tetrahedron ;

//----------------------------------

typedef Shape< 5 ,
        // Edges:
        TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 0 > ,
        TypeList< ShapeBoundary< Edge , 0 , 3 > ,
        TypeList< ShapeBoundary< Edge , 0 , 4 > ,
        TypeList< ShapeBoundary< Edge , 1 , 4 > ,
        TypeList< ShapeBoundary< Edge , 2 , 4 > ,
        TypeList< ShapeBoundary< Edge , 3 , 4 > ,
        TypeListEnd > > > > > > > > ,
        // Faces:
        TypeList< ShapeBoundary< TriangleFace , 0 , 1 , 4 > ,
        TypeList< ShapeBoundary< TriangleFace , 1 , 2 , 4 > ,
        TypeList< ShapeBoundary< TriangleFace , 2 , 3 , 4 > ,
        TypeList< ShapeBoundary< TriangleFace , 3 , 0 , 4 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 3 , 2 , 1 > ,
        TypeListEnd > > > > > ,
        3 , 3 > Pyramid ;

//----------------------------------

typedef Shape< 6 ,
        // Edges:
        TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 0 > ,
        TypeList< ShapeBoundary< Edge , 3 , 4 > ,
        TypeList< ShapeBoundary< Edge , 4 , 5 > ,
        TypeList< ShapeBoundary< Edge , 5 , 3 > ,
        TypeList< ShapeBoundary< Edge , 0 , 3 > ,
        TypeList< ShapeBoundary< Edge , 1 , 4 > ,
        TypeList< ShapeBoundary< Edge , 2 , 5 > ,
        TypeListEnd > > > > > > > > > ,
        // Faces:
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 1 , 4 , 3 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 1 , 2 , 5 , 4 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 3 , 5 , 2 > ,
        TypeList< ShapeBoundary< TriangleFace , 0 , 2 , 1 > ,
        TypeList< ShapeBoundary< TriangleFace , 3 , 4 , 5 > ,
        TypeListEnd > > > > > ,
        3 , 3 > Wedge ;

//----------------------------------

typedef Shape< 8 ,
        // Edges
        TypeList< ShapeBoundary< Edge , 0 , 1 > ,
        TypeList< ShapeBoundary< Edge , 1 , 2 > ,
        TypeList< ShapeBoundary< Edge , 2 , 3 > ,
        TypeList< ShapeBoundary< Edge , 3 , 0 > ,
        TypeList< ShapeBoundary< Edge , 4 , 5 > ,
        TypeList< ShapeBoundary< Edge , 5 , 6 > ,
        TypeList< ShapeBoundary< Edge , 6 , 7 > ,
        TypeList< ShapeBoundary< Edge , 7 , 4 > ,
        TypeList< ShapeBoundary< Edge , 0 , 4 > ,
        TypeList< ShapeBoundary< Edge , 1 , 5 > ,
        TypeList< ShapeBoundary< Edge , 2 , 6 > ,
        TypeList< ShapeBoundary< Edge , 3 , 7 > ,
        TypeListEnd > > > > > > > > > > > > ,
        // Faces:
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 1 , 5 , 4 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 1 , 2 , 6 , 5 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 2 , 3 , 7 , 6 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 4 , 7 , 3 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 3 , 2 , 1 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 4 , 5 , 6 , 7 > ,
        TypeListEnd > > > > > > ,
        3 , 3 > Hexahedron ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Special element shapes (bar and shell)

typedef Shape< 2, TypeList< ShapeBoundary< Edge , 0 , 1 > , TypeListEnd > ,
                  TypeListEnd ,
                  2 , 3 > LineBar ;

//----------------------------------

typedef TypeList< ShapeBoundary< Edge, 0, 1 >,
        TypeList< ShapeBoundary< Edge, 1, 0 >,
        TypeListEnd> > Line2DShellEdges ;

typedef Shape< 2, Line2DShellEdges , Line2DShellEdges ,  2 , 2 > Line2DShell ;

//----------------------------------

typedef Shape< 3, TriangleEdges,
                  TypeList< ShapeBoundary< TriangleFace , 0 , 1 , 2 > ,
                  TypeList< ShapeBoundary< TriangleFace , 0 , 2 , 1 > ,
                  TypeListEnd > > ,
                  3 , 3 > TriangleShell ;

//----------------------------------

typedef Shape< 4,
        QuadrilateralEdges,
        // Faces:
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 1 , 2 , 3 > ,
        TypeList< ShapeBoundary< QuadrilateralFace , 0 , 3 , 2 , 1 > ,
        TypeListEnd > > ,
        3 , 3 >  QuadrilateralShell ;

//----------------------------------------------------------------------

} // namespace element
} // namespace phdmesh

//----------------------------------------------------------------------

#endif

