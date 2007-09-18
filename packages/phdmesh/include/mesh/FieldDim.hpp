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

#ifndef phdmesh_FieldDim_hpp
#define phdmesh_FieldDim_hpp

//----------------------------------------------------------------------

#include <mesh/Types.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

class FieldDimension {
private:
  enum { MaxDim = MaximumFieldDimension };
  enum { OK = StaticAssert< MaxDim == 8 >::OK };
  enum { I_Length = MaxDim + 1 };
  enum { I_Size   = MaxDim + 2 };
  enum { I_NDim   = MaxDim + 3 };
  enum { MaxInfo  = MaxDim + 4 };
  unsigned m_info[ MaxInfo ];

public:

  FieldDimension();
  FieldDimension( const FieldDimension & rhs );
  FieldDimension & operator = ( const FieldDimension & rhs );

  FieldDimension( unsigned scalar_size ,
                  unsigned n0 ,     unsigned n1 = 0 , unsigned n2 = 0 ,
                  unsigned n3 = 0 , unsigned n4 = 0 , unsigned n5 = 0 ,
                  unsigned n6 = 0 , unsigned n7 = 0 );

  bool operator == ( const FieldDimension & rhs ) const ;
  bool operator != ( const FieldDimension & rhs ) const ;

  unsigned number_of_dimensions() const { return m_info[ I_NDim ]; }

  /** Individual dimension */
  unsigned operator[]( unsigned i ) const
    { return i < m_info[ I_NDim ] ? m_info[i+1] / m_info[i] : 0 ; }

  /** Length = product of dimensions */
  unsigned length() const { return m_info[ I_Length ]; }

  /** Size = byte size */
  unsigned size() const { return m_info[ I_Size ]; }

  unsigned offset( unsigned * const i )
    {
      unsigned n = 0 , j = m_info[ I_NDim ] ;
      while ( j ) { --j ; n += m_info[j] * i[j] ; } 
      return n ;
    }

  unsigned offset( unsigned i0 ) const
    { return i0 ; }

  unsigned offset( unsigned i0 , unsigned i1 ) const
    { return i0 + m_info[1] * i1 ; }

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 ; }

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ,
                   unsigned i3 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 + m_info[3] * i3 ; }

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ,
                   unsigned i3 , unsigned i4 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 + m_info[3] * i3 +
                  m_info[4] * i4 ; }

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ,
                   unsigned i3 , unsigned i4 , unsigned i5 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 + m_info[3] * i3 +
                  m_info[4] * i4 + m_info[5] * i5 ; }

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ,
                   unsigned i3 , unsigned i4 , unsigned i5 ,
                   unsigned i6 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 + m_info[3] * i3 +
                  m_info[4] * i4 + m_info[5] * i5 + m_info[6] * i6 ; };

  unsigned offset( unsigned i0 , unsigned i1 , unsigned i2 ,
                   unsigned i3 , unsigned i4 , unsigned i5 ,
                   unsigned i6 , unsigned i7 ) const
    { return i0 + m_info[1] * i1 + m_info[2] * i2 + m_info[3] * i3 +
                  m_info[4] * i4 + m_info[5] * i5 + m_info[6] * i6 +
                  m_info[7] * i7 ; }

  /** All dimensions, return number of dimensions */
  unsigned dimension( unsigned * const ) const ;

  /** Map offset to indices, return if offset is in range */
  bool indices( unsigned , unsigned * const ) const ;
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

