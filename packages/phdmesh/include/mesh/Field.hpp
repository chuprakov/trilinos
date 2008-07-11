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

#ifndef phdmesh_Field_hpp
#define phdmesh_Field_hpp

//----------------------------------------------------------------------

#include <iosfwd>
#include <string>
#include <vector>

#include <util/NumericEnum.hpp>
#include <util/FixedArray.hpp>
#include <util/CSet.hpp>

#include <mesh/Types.hpp>
#include <mesh/FieldTraits.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

/** Print the text name for a field, depending on the number of states.  */
std::ostream & operator << ( std::ostream & , const FieldBase & );

/** Print field and field dimension map entries on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const FieldBase & );

//----------------------------------------------------------------------
/** Field value states.
 *    Default state is 'None' or 'New' or 'N+1'
 */
enum FieldState {
  /* Exactly one state  */    StateNone = 0 ,
  /* Exactly two states */    StateNew  = 0 ,
                              StateOld  = 1 ,
  /* Three or more states */  StateNM1  = 1  /* N-1 */ ,
                              StateNM2  = 2  /* N-2 */ ,
                              StateNM3  = 3  /* N-3 */ ,
                              StateNM4  = 4  /* N-4 */ ,
                              StateNM5  = 5  /* N-5 */ };

enum { MaximumFieldStates = 6 };

const char * field_state_name( FieldState );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
class Field< void , DimensionTag ,
                    DimensionTag , DimensionTag ,
                    DimensionTag , DimensionTag ,
                    DimensionTag , DimensionTag > {
public:

  typedef void data_type ;

  Schema & schema() const { return m_schema ; }

  /** Volatile until schema is committed. */
  unsigned schema_ordinal() const { return m_schema_ordinal ; }

  const std::string & name() const { return m_name ; }

  template<typename NumType> bool type_is() const
    { return m_scalar_type == NumericEnum<NumType>::value ; }

  unsigned numeric_type_ordinal() const { return m_scalar_type ; }

  unsigned number_of_states() const { return m_num_states ; }

  FieldState state() const { return m_this_state ; }

  unsigned number_of_dimensions() const { return m_num_dim ; }

  const DimensionTag * const * dimension_tags() const
    { return m_dim_tags ; }

  unsigned max_size( EntityType ) const ;

  //----------------------------------------

  template<class A>
  CSet::Span<A> attribute() const { return m_cset.get<A>(); }

  //----------------------------------------
  /** An internal data structure that should never need to be
   *  used by a user of the phdMesh package.
   */
  struct Dim {
    const Part * part ;
    EntityType   type ;
    unsigned     stride[ MaximumFieldDimension ];
/*
    DimensionNoTag dim ;
*/

    Dim( EntityType t , const Part & p );
    Dim();
    Dim( const Dim & rhs );
    Dim & operator = ( const Dim & rhs );
  };

  /** Volatile until the schema is committed */
  const std::vector<Dim> & dimension() const ;

  /** Volatile until the schema is committed */
  const Dim & dimension( EntityType , const Part & ) const ;

  //----------------------------------------

private:

  friend class Schema ;

  template< typename Scalar , class Tag1 , class Tag2 ,
                              class Tag3 , class Tag4 ,
                              class Tag5 , class Tag6 ,
                              class Tag7 >
    friend class Field ;

  ~Field();

  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

  Field( Schema & ,
         const std::string & ,
         unsigned scalar_type ,
         const DimensionTag * ,
         const DimensionTag * ,
         const DimensionTag * ,
         const DimensionTag * ,
         const DimensionTag * ,
         const DimensionTag * ,
         const DimensionTag * ,
         unsigned number_of_states ,
         FieldState );

  std::vector<Dim> & dimension();

  CSet        m_cset ;
  std::string m_name ;
  Schema    & m_schema ;         // Mesh schema in which this field resides
  unsigned    m_schema_ordinal ; // Ordinal in the field set
  unsigned    m_scalar_type ;    // Ordinal in FieldTypes
  unsigned    m_num_dim ;        // Number of dimensions
  unsigned    m_num_states ;     // Number of states of this field
  FieldState  m_this_state ;     // Field state of this field

  std::vector<Dim>     m_dim_map ; // Only valid on StateNone
  Field              * m_field_states[ MaximumFieldStates ];
  const DimensionTag * m_dim_tags[ MaximumFieldDimension ];
};

//----------------------------------------------------------------------
// Maximum dimension is one less to allow introduction
// additional 'EntityDimension'.

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
                            class Tag5 , class Tag6 , class Tag7 >
class Field : public FieldBase {
public:

  typedef Scalar data_type ;

  typedef Tag1 dimension_tag_1 ;
  typedef Tag2 dimension_tag_2 ;
  typedef Tag3 dimension_tag_3 ;
  typedef Tag4 dimension_tag_4 ;
  typedef Tag5 dimension_tag_5 ;
  typedef Tag6 dimension_tag_6 ;
  typedef Tag7 dimension_tag_7 ;

  /** Dimension of field data for a single entity */
  typedef phdmesh::Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> Dimension ;

  /** Dimension of a contiguous block of field data for a
   *  homogeneous collection of entities. 
   */
  typedef typename
    DimensionAppend<Dimension,EntityDimension>::type BlockDimension ;

  /** Query field data dimension for a single entity. */
  Dimension dimension( EntityType entity_type , const Part & part ) const
    {
      const FieldBase::Dim & d = FieldBase::dimension( entity_type , part );
      return Dimension( d.stride );
    }

  /** Query this field for a given field state. */
  const Field & operator[]( FieldState state ) const
    { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

private:

  ~Field();
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

inline
FieldBase::Dim::Dim( EntityType t , const Part & p )
  : part( & p ) , type( t )
  { Copy<MaximumFieldDimension>( stride , 0u ); }

inline
FieldBase::Dim::Dim()
  : part( NULL ), type( EntityTypeEnd )
  { Copy<MaximumFieldDimension>( stride , 0u ); }

inline
FieldBase::Dim::Dim( const FieldBase::Dim & rhs )
  : part( rhs.part ), type( rhs.type )
  { Copy< MaximumFieldDimension >( stride , rhs.stride ); }

inline
FieldBase::Dim &
FieldBase::Dim::operator = ( const FieldBase::Dim & rhs )
  {
    part = rhs.part ;
    type = rhs.type ;
    Copy< MaximumFieldDimension >( stride , rhs.stride );
    return *this ;
  }

//----------------------------------------------------------------------
/** A field relation for fields that are a pointers to other fields.
 *
 *  If Entity 'e1' has a relation to Entity 'e2' that is in the
 *     domain of the relation stencil 'm_function' AND
 *     field 'm_root'   has a pointer scalar type 'T *' AND
 *     field 'm_target' has a scalar type 'T' AND
 *     field 'm_root'   has field data for Entity 'e1' AND
 *     field 'm_target' has field data for Entity 'e2' AND
 *     the 'e1' to 'e2' relation identifier is within the
 *     the 'm_root' field data size
 *   then
 *      field_data(*m_root,e1)[index] == field_data(*m_target,e2)
 *      where index = (*m_function)( e1.entity_type() ,
 *                                   e2.entity_type() ,
 *                                   relation_identifier_e1_to_e2 ,
 *                                   relation_kind_e1_to_e2 )
 *
 *  This data structure is used internally and should never need to be
 *  used by a user of the phdMesh package.
 */
struct FieldRelation {
  FieldBase          * m_root ;
  FieldBase          * m_target ;
  relation_stencil_ptr m_function ;

  FieldRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  FieldRelation( const FieldRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  FieldRelation & operator = ( const FieldRelation & rhs )
    {
      m_root = rhs.m_root ;
      m_target = rhs.m_target ;
      m_function = rhs.m_function ;
      return *this ;
    }
};

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

