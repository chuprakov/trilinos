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

#ifndef phdmesh_FieldData_hpp
#define phdmesh_FieldData_hpp

//----------------------------------------------------------------------

#include <mesh/Field.hpp>
#include <mesh/Kernel.hpp>
#include <mesh/Entity.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
/** Check validity of field and kernel: compatibility and existence.
 *  For performance none of the remaining field_data functions have
 *  internal validity checks.
 */
bool field_data_valid( const FieldBase & f ,
                       const Kernel & k ,
                       unsigned ord = 0,
                       const char * required_by = NULL );

inline
bool field_data_valid( const FieldBase & f ,
                       const Entity & e ,
                       const char * required_by = NULL )
{ return field_data_valid( f, e.kernel(), e.kernel_ordinal(), required_by ); }

//----------------------------------------------------------------------

template< class field_type >
inline
typename field_type::Dimension
field_dimension( const field_type & f , const Kernel & k )
{
  typedef typename field_type::Dimension Dim ;
  const Kernel::DataMap & pd = k.m_field_map[ f.schema_ordinal() ];
  return Dim( pd.m_stride , "phdmesh::field_dimension" );
}

template< class field_type >
inline
typename field_type::data_type *
field_data( const field_type & f , const Kernel & k )
{
  typedef unsigned char                  * byte_p ;
  typedef typename field_type::data_type * data_p ;

  data_p ptr = NULL ;

  const Kernel::DataMap & pd = k.m_field_map[ f.schema_ordinal() ];

  if ( pd.m_size ) {
    ptr = reinterpret_cast<data_p>(
          reinterpret_cast<byte_p>(k.m_entities) + pd.m_base );
  }
  return ptr ;
}

template< class field_type >
inline
typename field_type::data_type *
field_data( const field_type & f , const Kernel & k , unsigned i )
{
  typedef unsigned char                  * byte_p ;
  typedef typename field_type::data_type * data_p ;

  data_p ptr = NULL ;

  const Kernel::DataMap & pd = k.m_field_map[ f.schema_ordinal() ];

  if ( pd.m_size ) {
    ptr = reinterpret_cast<data_p>(
          reinterpret_cast<byte_p>(k.m_entities) + pd.m_base + pd.m_size * i );
  }
  return ptr ;
}

inline
unsigned field_data_size( const FieldBase & f , const Kernel & k )
{
  const Kernel::DataMap & pd = k.m_field_map[ f.schema_ordinal() ];
  return pd.m_size ;
}

//----------------------------------------------------------------------

template< class field_type >
inline
typename field_type::Dimension
field_dimension( const field_type & f , const Entity & e )
{ return field_dimension( f , e.kernel() ); }

template< class field_type >
inline
typename field_type::data_type *
field_data( const field_type & f , const Entity & e )
{ return field_data( f , e.kernel() , e.kernel_ordinal() ); }

inline
unsigned field_data_size( const FieldBase & f , const Entity & e )
{ return field_data_size( f , e.kernel() ); }

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

