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

#ifndef phdmesh_Schema_hpp
#define phdmesh_Schema_hpp

//----------------------------------------------------------------------

#include <util/Parallel.hpp>
#include <mesh/Types.hpp>
#include <mesh/Part.hpp>
#include <mesh/Field.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

/** Parallel Heterogeneous Dynamic Mesh.
 *  An dynamic unstructured mesh of mesh entities with
 *  subsets of parts partitioned into homogeneous kernels.
 */

class Schema {
public:

  //------------------------------------
  // Predefined parts:

  /** Universal: superset of all other parts. */
  Part & universal_part() const { return const_cast<Part&>(m_universal_part); }

  /** Uses: used by the local processor, subset of 'universal_part'.
   *  The aura is a subset of 'universal_part \ uses_part'.
   */
  Part & uses_part() const { return *m_uses_part ; }

  /** Owned: owned by the local processor, subset of 'uses_part'.
   *  Each mesh entity is owned by exactly one processor.
   */
  Part & owns_part()  const { return *m_owns_part ; }

  //------------------------------------
  /** Get an existing part of the given name and type.
   *  Return NULL if not present and required_by == NULL.
   *  If required and not present then throws an exception
   *  with the 'required_by' text.
   */
  Part * get_part( const std::string & ,
                   const char * required_by = NULL ) const ;

  /** Query all parts of the mesh.
   *  Once the mesh is commited the parts will know
   *  their ordinals within this vector.
   */
  const PartSet & get_parts() const { return m_universal_part.subsets(); }

  /** Declare a part of the given name.
   *  Redeclaration returns the previously declared part.
   */
  Part & declare_part( const std::string & );

  /** Declare a part that is defined as the
   *  intersection of the given part set.
   */
  Part & declare_part( const PartSet & );

  /** Declare a superset-subset relationship */
  void declare_part_subset( Part & superset , Part & subset );

  /** Declare an attribute on a part */
  template<class T>
  CSet::Span<T> declare_part_attribute( Part & , const T * , bool );

  //------------------------------------
  /** Get a field, return NULL if it does not exist.
   *  If T != void or NDim != 0 then an exception will be thrown
   *  if the field exits and the type or number of dimensions does not match.
   *  If required and not present then throws an exception
   *  with the 'required_by' text.
   */
  template<typename T,unsigned NDim>
  Field<T,NDim> * get_field( EntityType entity_type ,
                             const std::string & name ,
                             const char * required_by = NULL ) const ;

  /** Get all fields associated with the given entity type */
  const std::vector< Field<void,0> *> & get_fields( EntityType t ) const
    { return m_fields[ t ]; }

  /** Declare a field within the mesh.
   *  Redeclaration with compatible parameters returns the
   *  previously declared field.
   *  Redeclaration with incompatible parameters throws an exception.
   */
  template<typename T,unsigned NDim>
  Field<T,NDim> & declare_field( EntityType entity_type ,
                                 const std::string & name ,
                                 unsigned number_of_states = 1 );

  /** Declare a field to have a dimension over a given part */
  void declare_field_dimension( Field<void,0> & field , const Part & ,
                                unsigned n0 ,     unsigned n1 = 0 ,
                                unsigned n2 = 0 , unsigned n3 = 0 ,
                                unsigned n4 = 0 , unsigned n5 = 0 ,
                                unsigned n6 = 0 , unsigned n7 = 0 );

  /** Declare an attribute on a field */
  template<class T>
  CSet::Span<T> declare_field_attribute( Field<void,0> & , const T * , bool );
  //------------------------------------
  /** Commit the part and field declarations.
   *  Verifies consistency and assigns ordinals for faster usage.
   *  No more declarations (parts, part-subsets, fields, field-dimensions)
   *  can be made.
   */
  void commit();

  //------------------------------------

  void assert_committed( const char * ) const ;

  void assert_not_committed( const char * ) const ;

  void assert_same_schema( const char * , const Schema & ) const ;

  bool is_commit() const { return m_commit ; }

  ~Schema();

  Schema();

private:

  Schema( const Schema & );
  Schema & operator = ( const Schema & );

  bool   m_commit ;
  Part   m_universal_part ;
  Part * m_uses_part ;
  Part * m_owns_part ;

  std::vector< Field<void,0> * > m_fields[ EntityTypeMaximum ];

  Field<void,0> & declare_field( EntityType ,
                                 const std::string & ,
                                 unsigned arg_scalar_type ,
                                 unsigned arg_num_dim ,
                                 unsigned arg_num_states );

  Field<void,0> * get_field( bool ,
                             EntityType ,
                             const std::string & ,
                             unsigned arg_scalar_type ,
                             unsigned arg_num_dim ,
                             unsigned arg_num_states ,
                             const char * required_by ) const ;

  void clean_field_dimension();
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

template<class T>
inline
CSet::Span<T>
Schema::declare_part_attribute( Part & p , const T * a , bool d )
{
  assert_not_committed( "phdmesh::Schema::declare_part_attribute" );
  return p.m_cset.insert<T>( a , d );
}

template<class T>
inline
CSet::Span<T>
Schema::declare_field_attribute( Field<void,0> & f , const T * a , bool d )
{
  assert_not_committed( "phdmesh::Schema::declare_field_attribute" );
  return f.m_cset.insert<T>( a , d );
}

template<typename T,unsigned NDim>
inline
Field<T,NDim> *
Schema::get_field( EntityType entity_type ,
                   const std::string & name ,
                   const char * required_by ) const
{
  enum { ftype = NumericEnum<T>::value };

  Field<void,0> * const f =
    get_field( false , entity_type , name , ftype , NDim , 0 , required_by );

  return static_cast< Field<T,NDim> *>( f );
}

template<typename T,unsigned NDim>
inline
Field<T,NDim> &
Schema::declare_field( EntityType entity_type ,
                       const std::string & name ,
                       unsigned number_of_states )
{
  enum { ftype = NumericEnum<T>::value };

  Field<void,0> & f =
    declare_field( entity_type, name, ftype, NDim, number_of_states );

  return static_cast<Field<T,NDim> &>( f );
}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

