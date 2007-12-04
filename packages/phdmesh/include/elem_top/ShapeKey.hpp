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

#ifndef phdmesh_ElementShapeKey_hpp
#define phdmesh_ElementShapeKey_hpp

namespace phdmesh {
namespace element {

//----------------------------------------------------------------------

struct ShapeKey {

  template< class Shape >
  struct Encoding {
    enum { value = ( ((unsigned) Shape::number_vertex ) << 24 ) |
                   ( ((unsigned) Shape::number_edge   ) << 16 ) |
                   ( ((unsigned) Shape::number_side   ) <<  8 ) |
                   ( ((unsigned) Shape::boundary      ) <<  0 ) };
  };

  unsigned char m_number_vertex ;
  unsigned char m_number_edge ;
  unsigned char m_number_side ;
  unsigned char m_boundary ;

  ~ShapeKey() {}

  ShapeKey()
    : m_number_vertex( 0 ),
      m_number_edge( 0 ),
      m_number_side( 0 ),
      m_boundary( 0 ) {}

  ShapeKey( const ShapeKey & rhs )
    : m_number_vertex( rhs.m_number_vertex ),
      m_number_edge( rhs.m_number_edge ),
      m_number_side( rhs.m_number_side ),
      m_boundary( rhs.m_boundary ) {}

  ShapeKey( unsigned char arg_number_vertex ,
       unsigned char arg_number_edge ,
       unsigned char arg_number_side ,
       unsigned char arg_boundary )
    : m_number_vertex( arg_number_vertex ),
      m_number_edge( arg_number_edge ),
      m_number_side( arg_number_side ),
      m_boundary( arg_boundary ) {}

  explicit
  ShapeKey( unsigned encoded )
    : m_number_vertex( ( encoded >> 24 ) & 0x0f ),
      m_number_edge(   ( encoded >> 16 ) & 0x0f ),
      m_number_side(   ( encoded >>  8 ) & 0x0f ),
      m_boundary(      ( encoded >>  0 ) & 0x0f ) {}

  unsigned encoding() const
    {
      return ( ((unsigned) m_number_vertex ) << 24 ) |
             ( ((unsigned) m_number_edge   ) << 16 ) |
             ( ((unsigned) m_number_side   ) <<  8 ) |
             ( ((unsigned) m_boundary      ) <<  0 ) ;
    }

  ShapeKey & operator = ( const ShapeKey & rhs )
    {
      m_number_vertex = rhs.m_number_vertex ;
      m_number_edge   = rhs.m_number_edge ;
      m_number_side   = rhs.m_number_side ;
      m_boundary      = rhs.m_boundary ;
      return *this ;
    }
};

} // namespace element
} // namespace phdmesh

//----------------------------------------------------------------------

#endif

