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

#ifndef phdmesh_FieldDeclarations_hpp
#define phdmesh_FieldDeclarations_hpp

#include <mesh/FieldTraits.hpp>
#include <mesh/MetaData.hpp>

namespace phdmesh {

typedef Field<double>                  ScalarField ;
typedef Field<double, Cartesian>       VectorField ;
typedef Field<double, FullTensor>      FullTensorField ;
typedef Field<double, SymmetricTensor> SymmetricTensorField ;

//----------------------------------------------------------------------
// Declaring fields of these fundemental types:

ScalarField &
declare_scalar_field_on_all_nodes( MeshMetaData & , const std::string & );

ScalarField &
declare_scalar_field_on_all_elements( MeshMetaData & , const std::string & );

VectorField &
declare_vector_field_on_all_nodes( MeshMetaData & , const std::string & , unsigned );

VectorField &
declare_vector_field_on_all_elements( MeshMetaData & , const std::string & , unsigned );

FullTensorField &
declare_full_tensor_field_on_all_nodes( MeshMetaData & , const std::string & , unsigned );

FullTensorField &
declare_full_tensor_field_on_all_elements( MeshMetaData & , const std::string & , unsigned );

SymmetricTensorField &
declare_symmetric_tensor_field_on_all_nodes( MeshMetaData & , const std::string & , unsigned );

SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements( MeshMetaData & , const std::string & , unsigned );

//----------------------------------------------------------------------
// Scalar fields:

template< typename Type >
Field<Type> & put_field_on_all_nodes( Field<Type> & );

template< typename Type >
Field<Type> & put_field_on_all_elements( Field<Type> & );

template< typename Type >
Field<Type> & put_field_on_nodes( Field<Type> & , Part & );

template< typename Type >
Field<Type> & put_field_on_elements( Field<Type> &, Part & );

//----------------------------------------------------------------------
// 1D Array fields:

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_all_nodes( Field<Type,T1> & , unsigned n1 );

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_all_elements( Field<Type,T1> & , unsigned n1 );

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_nodes( Field<Type,T1> & , Part & , unsigned n1 );

template< typename Type , class T1 >
Field<Type,T1> & put_field_on_elements( Field<Type,T1> &, Part &, unsigned n1 );

//----------------------------------------------------------------------
// 2D Array fields:

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_all_nodes( Field<Type,T1,T2> &, unsigned n1 , unsigned n2 );

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_all_elements( Field<Type,T1,T2> &, unsigned n1, unsigned n2);

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_nodes( Field<Type,T1,T2> &, Part &, unsigned n1, unsigned n2);

template< typename Type , class T1 , class T2 >
Field<Type,T1,T2> &
put_field_on_elements( Field<Type,T1,T2> &, Part &, unsigned n1, unsigned n2 );

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Template implementations follow:
//----------------------------------------------------------------------

namespace phdmesh {

template< typename Type >
inline
Field<Type> & put_field_on_nodes( Field<Type> & f , Part & p )
{ f.mesh_meta_data().put_field( f , Node , p );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> &
put_field_on_nodes( Field<Type,T1> & f , Part & p , unsigned n1 )
{ f.mesh_meta_data().put_field( f , Node , p , n1 );  return f ; }

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> &
put_field_on_nodes( Field<Type,T1,T2> & f , Part & p ,
                    unsigned n1 , unsigned n2 )
{ f.mesh_meta_data().put_field( f , Node , p , n1 , n2 );  return f ; }

//----------------------------------------------------------------------

template< typename Type >
inline
Field<Type> & put_field_on_elements( Field<Type> & f , Part & p )
{ f.mesh_meta_data().put_field( f , Element , p );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> &
put_field_on_elements( Field<Type,T1> & f , Part & p , unsigned n1 )
{ f.mesh_meta_data().put_field( f , Element , p , n1 );  return f ; }

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> &
put_field_on_elements( Field<Type,T1,T2> & f , Part & p ,
                                unsigned n1 , unsigned n2 )
{ f.mesh_meta_data().put_field( f , Element , p , n1 , n2 );  return f ; }

//----------------------------------------------------------------------

template< typename Type >
inline
Field<Type> & put_field_on_all_nodes( Field<Type> & f )
{ put_field_on_nodes( f , f.mesh_meta_data().universal_part() );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> & put_field_on_all_nodes( Field<Type,T1> & f , unsigned n1 )
{
  put_field_on_nodes( f , f.mesh_meta_data().universal_part() , n1 );
  return f ;
}

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> &
put_field_on_all_nodes( Field<Type,T1,T2> & f , unsigned n1 , unsigned n2 )
{
  put_field_on_nodes( f , f.mesh_meta_data().universal_part() , n1 , n2 );
 return f ;
}

//----------------------------------------------------------------------

template< typename Type >
inline
Field<Type> & put_field_on_all_elements( Field<Type> & f )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() );  return f ; }

template< typename Type , class T1 >
inline
Field<Type,T1> & put_field_on_all_elements( Field<Type,T1> & f , unsigned n1 )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() , n1 );  return f ; }

template< typename Type , class T1 , class T2 >
inline
Field<Type,T1,T2> & put_field_on_all_elements( Field<Type,T1,T2> & f ,
                                unsigned n1 , unsigned n2 )
{ put_field_on_elements( f , f.mesh_meta_data().universal_part() , n1 , n2 );  return f ; }

//----------------------------------------------------------------------

inline
ScalarField &
declare_scalar_field_on_all_nodes( MeshMetaData & md , const std::string & n )
{
  return md.put_field( md.declare_field<ScalarField>(n) ,
                       Node , md.universal_part() );
}

inline
ScalarField &
declare_scalar_field_on_all_elements( MeshMetaData & md ,
                                      const std::string & n )
{
  return md.put_field( md.declare_field<ScalarField>(n) ,
                       Element , md.universal_part() );
}

inline
VectorField &
declare_vector_field_on_all_nodes(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<VectorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
VectorField &
declare_vector_field_on_all_elements(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<VectorField>(s),
                       Element , md.universal_part() , n1 );
}

inline
FullTensorField &
declare_full_tensor_field_on_all_nodes(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<FullTensorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
FullTensorField &
declare_full_tensor_field_on_all_elements(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<FullTensorField>(s),
                       Element , md.universal_part() , n1 );
}

inline
SymmetricTensorField &
declare_symmetric_tensor_field_on_all_nodes(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<SymmetricTensorField>(s),
                       Node , md.universal_part() , n1 );
}

inline
SymmetricTensorField &
declare_symmetric_tensor_field_on_all_elements(
  MeshMetaData & md , const std::string & s , unsigned n1 )
{
  return md.put_field( md.declare_field<SymmetricTensorField>(s) ,
                       Element , md.universal_part() , n1 );
}

}

#endif

