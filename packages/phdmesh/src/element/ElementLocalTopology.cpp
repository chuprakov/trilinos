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

#include <stdlib.h>
#include <element/LocalTopology.hpp>
#include <element/TriangleTopology.hpp>
#include <element/QuadrilateralTopology.hpp>
#include <element/TetrahedronTopology.hpp>
#include <element/PyramidTopology.hpp>
#include <element/WedgeTopology.hpp>
#include <element/HexahedronTopology.hpp>

namespace phdmesh {

namespace {

template< class TopTraits > const LocalTopology * top_descriptor();
template< class IList > struct IndexListData ;

//----------------------------------------------------------------------

template<>
const LocalTopology * top_descriptor<TypeListEnd>()
{ return (const LocalTopology *) NULL ; }

template<>
struct IndexListData< TypeListEnd > { static const unsigned * data(); };

const unsigned * IndexListData< TypeListEnd >::data()
{ return (const unsigned *) NULL ; }

//----------------------------------------------------------------------

template< class TopTraits >
const LocalTopology * top_descriptor()
{ return TopTraits::descriptor(); }

template< unsigned  I0 , unsigned  I1 , unsigned  I2 , unsigned  I3 ,
          unsigned  I4 , unsigned  I5 , unsigned  I6 , unsigned  I7 ,
          unsigned  I8 , unsigned  I9 , unsigned I10 , unsigned I11 ,
          unsigned I12 , unsigned I13 , unsigned I14 , unsigned I15 ,
          unsigned I16 , unsigned I17 , unsigned I18 , unsigned I19 ,
          unsigned I20 , unsigned I21 , unsigned I22 , unsigned I23 ,
          unsigned I24 , unsigned I25 , unsigned I26 , unsigned I27 ,
          unsigned I28 , unsigned I29 , unsigned I30 , unsigned I31 >
struct IndexListData< IndexList<  I0 ,  I1 ,  I2 ,  I3 ,
                                  I4 ,  I5 ,  I6 ,  I7 ,
                                  I8 ,  I9 , I10 , I11 ,
                                 I12 , I13 , I14 , I15 ,
                                 I16 , I17 , I18 , I19 ,
                                 I20 , I21 , I22 , I23 ,
                                 I24 , I25 , I26 , I27 ,
                                 I28 , I29 , I30 , I31 > >
{
  static const unsigned * data()
  {
    static const unsigned self[] = { I0 ,  I1 ,  I2 ,  I3 ,
                                     I4 ,  I5 ,  I6 ,  I7 ,
                                     I8 ,  I9 , I10 , I11 ,
                                    I12 , I13 , I14 , I15 ,
                                    I16 , I17 , I18 , I19 ,
                                    I20 , I21 , I22 , I21 ,
                                    I24 , I25 , I26 , I27 ,
                                    I28 , I29 , I30 , I31 };
    return self ;
  }
};

//----------------------------------------------------------------------

#define BOUNDARY_DATA_ENTRY( N )	\
  {	\
    top_descriptor< typename TypeListAt< BTypeList , N >::type >() ,	\
    IndexListData<  typename TypeListAt< BNodeMap ,  N >::type >::data() , \
    IndexListData<  typename TypeListAt< BNodeMap ,  N >::type >::data() , \
  }

template< class BTypeList , class BNodeMap >
const LocalTopology::Boundary * boundary_data()
{
  enum { Maximum = 16 , Requested = TypeListLength< BTypeList >::value };

  StaticAssert< Requested <= Maximum >::ok();

  static LocalTopology::Boundary data[ Maximum ] = {
    BOUNDARY_DATA_ENTRY(  0 ) ,
    BOUNDARY_DATA_ENTRY(  1 ) ,
    BOUNDARY_DATA_ENTRY(  2 ) ,
    BOUNDARY_DATA_ENTRY(  3 ) ,
    BOUNDARY_DATA_ENTRY(  4 ) ,
    BOUNDARY_DATA_ENTRY(  5 ) ,
    BOUNDARY_DATA_ENTRY(  6 ) ,
    BOUNDARY_DATA_ENTRY(  7 ) ,
    BOUNDARY_DATA_ENTRY(  8 ) ,
    BOUNDARY_DATA_ENTRY(  9 ) ,
    BOUNDARY_DATA_ENTRY( 10 ) ,
    BOUNDARY_DATA_ENTRY( 11 ) ,
    BOUNDARY_DATA_ENTRY( 12 ) ,
    BOUNDARY_DATA_ENTRY( 13 ) ,
    BOUNDARY_DATA_ENTRY( 14 ) ,
    BOUNDARY_DATA_ENTRY( 15 )
  };

  return data ;
}

#undef BOUNDARY_DATA_ENTRY

//----------------------------------------------------------------------

template< class Traits >
const LocalTopology * local_topology_descriptor()
{
  static LocalTopology self = {
    Traits::topological_rank ,
    Traits::minimum_dimension ,
    Traits::maximum_dimension ,
    Traits::number_vertex ,
    Traits::number_edge ,
    Traits::number_side ,
    Traits::number_node ,
    Traits::key ,
    Traits::is_boundary ,
    boundary_data< typename Traits::BoundaryEdgeTypeList ,
                   typename Traits::BoundaryEdgeNodeMap >() ,
    boundary_data< typename Traits::BoundarySideTypeList ,
                   typename Traits::BoundarySideNodeMap >() };

  return & self ;
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#define DECLARE_TOPOLOGY_DESCRIPTOR( T )	\
  const LocalTopology * T::descriptor()	\
  { return local_topology_descriptor< T >(); }

DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryPoint<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryPoint<1> )

DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryEdge<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryEdge<2> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryEdge<3> )

DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryTriangle<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryTriangle<3> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryTriangle<6> )

DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryQuadrilateral<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryQuadrilateral<4> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryQuadrilateral<8> )
DECLARE_TOPOLOGY_DESCRIPTOR( BoundaryQuadrilateral<9> )

//----------------------------------------------------------------------

DECLARE_TOPOLOGY_DESCRIPTOR( Sphere<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Sphere<1> )

//----------------------------------------------------------------------

DECLARE_TOPOLOGY_DESCRIPTOR( Line<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Line<2> )
DECLARE_TOPOLOGY_DESCRIPTOR( Line<3> )

//----------------------------------------------------------------------

DECLARE_TOPOLOGY_DESCRIPTOR( Triangle<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Triangle<3> )
DECLARE_TOPOLOGY_DESCRIPTOR( Triangle<6> )

DECLARE_TOPOLOGY_DESCRIPTOR( Quadrilateral<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Quadrilateral<4> )
DECLARE_TOPOLOGY_DESCRIPTOR( Quadrilateral<8> )
DECLARE_TOPOLOGY_DESCRIPTOR( Quadrilateral<9> )

//----------------------------------------------------------------------

DECLARE_TOPOLOGY_DESCRIPTOR( Tetrahedron<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Tetrahedron<4> )

DECLARE_TOPOLOGY_DESCRIPTOR( Pyramid<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Pyramid<5> )

DECLARE_TOPOLOGY_DESCRIPTOR( Wedge<0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Wedge<6> )

DECLARE_TOPOLOGY_DESCRIPTOR( Hexahedron< 0> )
DECLARE_TOPOLOGY_DESCRIPTOR( Hexahedron< 8> )
DECLARE_TOPOLOGY_DESCRIPTOR( Hexahedron<20> )
DECLARE_TOPOLOGY_DESCRIPTOR( Hexahedron<27> )

//----------------------------------------------------------------------

}

