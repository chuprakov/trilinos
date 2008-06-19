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

#ifndef phdmesh_LocalTopology_hpp
#define phdmesh_LocalTopology_hpp

#include <util/TypeList.hpp>
#include <util/IndexList.hpp>

namespace phdmesh {

// Element boundary local topologies:

template< unsigned NumNode = 0 > struct BoundaryPoint ;
template< unsigned NumNode = 0 > struct BoundaryEdge ;

// Element local topologies:

template< unsigned NumNode = 0 > struct Sphere ;
template< unsigned NumNode = 0 > struct Line ;         // 1D only

//----------------------------------------------------------------------
/** Traits for an element local topology.
 *
 *  Top::Shape                        The topology traits without nodes
 *  Top::key                          Unique enumeration key for the topology
 *  Top::number_vertex
 *  Top::number_edge
 *  Top::number_side                  duplicates edges in 2D
 *  Top::number_node
 *  Top::is_boundary                  if describing an element boundary
 *  Top::topological_rank
 *  Top::mimimum_dimension            minimum valid spatial dimension
 *  Top::maximum_dimension            maximum valid spatial dimension
 * 
 *  Top::edge<I>                      defined for I < number_edge
 *  Top::edge<I>::traits              local topology traits for edge #I
 *  Top::edge<I>::vertex<J>::ordinal  if J < edge<I>::traits::number_vertex
 *  Top::edge<I>::node<J>::ordinal    if J < edge<I>::traits::number_node
 * 
 *  Top::side<I>                      defined for I < number_side
 *  Top::side<I>::traits              local topology traits for side #I
 *  Top::side<I>::vertex<J>::ordinal  if J < side<I>::traits::number_vertex
 *  Top::side<I>::node<J>::ordinal    if J < side<I>::traits::number_node
 *
 *  const LocalTopology * Top::descriptor()       
 *
 *  The first 'number_vertex' nodes must be placed at the vertices.
 */
template< unsigned Topological_Rank ,
          unsigned Minimum_Dimension ,
          unsigned Maximum_Dimension ,
          unsigned Number_Vertex ,
          unsigned Number_Node  = 0 ,
          class    EdgeTypeList = TypeListEnd ,
          class    EdgeNodeMap  = TypeListEnd ,
          class    SideTypeList = TypeListEnd ,
          class    SideNodeMap  = TypeListEnd >
struct LocalTopologyTraits {
private:
  enum { edge_types_length = TypeListLength< EdgeTypeList >::value ,
         edge_maps_length  = TypeListLength< EdgeNodeMap  >::value ,
         side_types_length = TypeListLength< SideTypeList >::value ,
         side_maps_length  = TypeListLength< SideNodeMap  >::value };

  enum { OK_edge = StaticAssert< edge_types_length == edge_maps_length >:: OK };

  enum { OK_side = StaticAssert< side_types_length == side_maps_length >:: OK };

  enum { OK_key  = StaticAssert< ( 0 == Topological_Rank  >> 4 ) &&
                                 ( 0 == Number_Vertex     >> 4 ) &&
                                 ( 0 == edge_types_length >> 4 ) &&
                                 ( 0 == side_types_length >> 4 ) &&
                                 ( 0 == Number_Node       >> 16 ) >::OK };
public:

  enum { topological_rank  = Topological_Rank ,
         minimum_dimension = Minimum_Dimension ,
         maximum_dimension = Maximum_Dimension ,
         number_vertex     = Number_Vertex ,
         number_edge       = edge_types_length ,
         number_side       = side_types_length ,
         number_node       = Number_Node };

  enum { key = ( Topological_Rank  << 28 ) |
               ( Number_Vertex     << 24 ) |
               ( edge_types_length << 20 ) |
               ( side_types_length << 16 ) |
               ( Number_Node ) };

  enum { is_boundary = Topological_Rank < Minimum_Dimension };

  //----------------------------------

  template< class BTypes , class BMaps , unsigned I >
  struct Boundary {
  private:
    enum { OK  = StaticAssert< I < TypeListLength<BMaps>::value >::OK };
    typedef typename TypeListAt< BMaps  , I >::type boundary_map ;
  public:

    typedef typename TypeListAt< BTypes , I >::type traits ;

    template< unsigned J >
    class vertex {
    private:
      enum { OK = StaticAssert< J < traits::number_vertex >::OK };
    public:
      enum { ordinal = IndexListAt< boundary_map , J >::value };
    };

    template< unsigned J >
    class node {
    private:
      enum { OK = StaticAssert< J < traits::number_node >::OK };
    public:
      enum { ordinal = IndexListAt< boundary_map , J >::value };
    };
  };

  //----------------------------------

  template< unsigned I >
  struct edge : public Boundary< EdgeTypeList , EdgeNodeMap , I > {};

  template< unsigned I >
  struct side : public Boundary< SideTypeList , SideNodeMap , I > {};

  //----------------------------------

  typedef EdgeTypeList BoundaryEdgeTypeList ;
  typedef SideTypeList BoundarySideTypeList ;
  typedef EdgeNodeMap  BoundaryEdgeNodeMap ;
  typedef SideNodeMap  BoundarySideNodeMap ;
};

//----------------------------------------------------------------------
/** The runtime version of the local topology information */

struct LocalTopology {
  const char * name ;
  unsigned topological_rank ;
  unsigned minimum_dimension ;
  unsigned maximum_dimension ;
  unsigned number_vertex ;
  unsigned number_edge ;
  unsigned number_side ;
  unsigned number_node ;
  unsigned key ;
  bool     is_boundary ;

  struct Boundary {
    const LocalTopology * traits ;
    const unsigned      * vertex ;
    const unsigned      * node ;
  };

  const Boundary * edge ;
  const Boundary * side ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
struct BoundaryPoint<0>
  : public LocalTopologyTraits< 0 , 1 , 3, // Dimension
                                0 , 0 > // #Vertex, #Node
{
  typedef BoundaryPoint<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryPoint<1>
  : public LocalTopologyTraits< 0 , 1 , 3, // Dimension
                                0 , 1 > // #Vertex, #Node
{
  typedef BoundaryPoint<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

template<>
struct BoundaryEdge<0>
  : public LocalTopologyTraits< 1 , 2 , 3, // Dimension
                                2 , 0 > // #Vertex, #Node
{
  typedef BoundaryEdge<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryEdge<2>
  : public LocalTopologyTraits< 1 , 2 , 3, // Dimension
                                2 , 2 > // #Vertex, #Node
{
  typedef BoundaryEdge<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct BoundaryEdge<3>
  : public LocalTopologyTraits< 1 , 2 , 3, // Dimension
                                2 , 3 > // #Vertex, #Node
{
  typedef BoundaryEdge<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------
// Solid element shapes:

template<>
struct Sphere<0>
  : public LocalTopologyTraits< 3 , 1 , 3, // Dimension
                                0 , 0 > // #Vertex, #Node
{
  typedef Sphere<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Sphere<1>
  : public LocalTopologyTraits< 3 , 1 , 3, // Dimension
                                0 , 1 > // #Vertex, #Node
{
  typedef Sphere<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------

template<>
struct Line<0>
  : public LocalTopologyTraits< 1 , 1 , 1, // Dimension
                                2 , 0 > // #Vertex, #Node
{
  typedef Line<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Line<2>
  : public LocalTopologyTraits< 1 , 1 , 1, // Dimension
                                2 , 2 > // #Vertex, #Node
{
  typedef Line<0> Shape ;
  static const LocalTopology * descriptor();
};

template<>
struct Line<3>
  : public LocalTopologyTraits< 1 , 1 , 1, // Dimension
                                2 , 3 > // #Vertex, #Node
{
  typedef Line<0> Shape ;
  static const LocalTopology * descriptor();
};

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

