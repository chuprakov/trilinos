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
   *  An exception will be thrown
   *  if the field exits and the type or number of dimensions does not match.
   *  If required and not present then throws an exception
   *  with the 'required_by' text.
   */
  template< class field_type >
  field_type * get_field( const std::string & name ,
                          const char * required_by = NULL ) const
    {
      typedef typename field_type::data_type Scalar ;
      typedef typename field_type::dimension_traits_1 Traits1 ;
      typedef typename field_type::dimension_traits_2 Traits2 ;
      typedef typename field_type::dimension_traits_3 Traits3 ;
      typedef typename field_type::dimension_traits_4 Traits4 ;
      typedef typename field_type::dimension_traits_5 Traits5 ;
      typedef typename field_type::dimension_traits_6 Traits6 ;
      typedef typename field_type::dimension_traits_7 Traits7 ;

      return static_cast< field_type * >(
        get_field_base( name ,
                        NumericEnum<Scalar>::value ,
                        Traits1::descriptor() ,
                        Traits2::descriptor() ,
                        Traits3::descriptor() ,
                        Traits4::descriptor() ,
                        Traits5::descriptor() ,
                        Traits6::descriptor() ,
                        Traits7::descriptor() ,
                        -1 , required_by ) );
    }

  /** Get all fields associated with the given entity type */
  const std::vector< FieldBase * > & get_fields() const
    { return m_fields ; }

  /** Declare a field within the mesh.
   *  Redeclaration with compatible parameters returns the
   *  previously declared field.
   *  Redeclaration with incompatible parameters throws an exception.
   */
  template< class field_type >
  field_type & declare_field( const std::string & name ,
                              unsigned number_of_states = 1 )
    {
      typedef typename field_type::data_type Scalar ;
      typedef typename field_type::dimension_traits_1 Traits1 ;
      typedef typename field_type::dimension_traits_2 Traits2 ;
      typedef typename field_type::dimension_traits_3 Traits3 ;
      typedef typename field_type::dimension_traits_4 Traits4 ;
      typedef typename field_type::dimension_traits_5 Traits5 ;
      typedef typename field_type::dimension_traits_6 Traits6 ;
      typedef typename field_type::dimension_traits_7 Traits7 ;

      return static_cast< field_type & >(
        declare_field_base( name ,
                            NumericEnum<Scalar>::value ,
                            Traits1::descriptor() ,
                            Traits2::descriptor() ,
                            Traits3::descriptor() ,
                            Traits4::descriptor() ,
                            Traits5::descriptor() ,
                            Traits6::descriptor() ,
                            Traits7::descriptor() ,
                            number_of_states ) );
    }

  /** Declare a field to have a size over a given entity type and part.
   */
  template< class field_type >
  void declare_field_size( field_type & arg_field ,
                           unsigned     arg_entity_type ,
                           const Part & arg_part ,
                           const typename field_type::Dimension & arg_dim )
    {
      declare_field_stride( arg_field , arg_entity_type , arg_part ,
                            arg_dim.stride );
    }

  template< class field_type >
  void declare_field_exists( field_type & arg_field ,
                             unsigned     arg_entity_type ,
                             const Part & arg_part )
    {
      declare_field_stride( arg_field , arg_entity_type , arg_part , NULL );
    }

  /** Declare an attribute on a field */
  template<class T>
  CSet::Span<T> declare_field_attribute( FieldBase & , const T * , bool );

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
  Part   m_universal_part ; /* Subset list contains all other parts */
  Part * m_uses_part ;
  Part * m_owns_part ;

  std::vector< FieldBase * > m_fields ;

  void declare_field_stride( FieldBase & ,
                             unsigned , const Part & ,
                             const unsigned * );
  
  FieldBase & declare_field_base( const std::string & ,
                                  unsigned arg_scalar_type ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  const DimensionTraits * ,
                                  unsigned arg_num_states );

  FieldBase * get_field_base( const std::string & ,
                              unsigned arg_scalar_type ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              const DimensionTraits * ,
                              int arg_num_states ,
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
  return p.m_cset.template insert<T>( a , d );
}

template<class T>
inline
CSet::Span<T>
Schema::declare_field_attribute( FieldBase & f , const T * a , bool d )
{
  assert_not_committed( "phdmesh::Schema::declare_field_attribute" );
  return f.m_cset.template insert<T>( a , d );
}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

