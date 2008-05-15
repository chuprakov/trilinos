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
class Field< void , DimensionTraits ,
                    DimensionTraits , DimensionTraits ,
                    DimensionTraits , DimensionTraits ,
                    DimensionTraits , DimensionTraits > {
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

  const DimensionTraits * const * dimension_traits() const
    { return m_dim_traits ; }

  unsigned max_size( unsigned entity_type ) const ;

  //----------------------------------------

  struct Dim {
    const Part * part ;
    unsigned     rank ;
    unsigned     stride[ MaximumFieldDimension ];

    Dim( unsigned r , const Part & p ) : part( & p ) , rank( r )
      { Copy<MaximumFieldDimension>( stride , 0u ); }

    Dim() : part( NULL ), rank( 0 )
      { Copy<MaximumFieldDimension>( stride , 0u ); }

    Dim( const Dim & rhs )
      : part( rhs.part ), rank( rhs.rank )
      { Copy< MaximumFieldDimension >( stride , rhs.stride ); }

    Dim & operator = ( const Dim & rhs )
      {
        part = rhs.part ;
        rank = rhs.rank ;
        Copy< MaximumFieldDimension >( stride , rhs.stride );
        return *this ;
      }
  };

  /** Volatile until the schema is committed */
  const std::vector<Dim> & dimension() const ;

  /** Volatile until the schema is committed */
  const Dim & dimension( unsigned , const Part & ) const ;

  //----------------------------------------

  template<class A>
  CSet::Span<A> attribute() const { return m_cset.get<A>(); }

  //----------------------------------------

private:

  friend class Schema ;

  template< typename Scalar ,
            class Traits1 , class Traits2 ,
            class Traits3 , class Traits4 ,
            class Traits5 , class Traits6 ,
            class Traits7 >
    friend class Field ;

  ~Field();

  Field();
  Field( const Field & );
  Field & operator = ( const Field & );

  Field( Schema & ,
         const std::string & ,
         unsigned scalar_type ,
         const DimensionTraits * ,
         const DimensionTraits * ,
         const DimensionTraits * ,
         const DimensionTraits * ,
         const DimensionTraits * ,
         const DimensionTraits * ,
         const DimensionTraits * ,
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

  std::vector<Dim>        m_dim_map ; // Only valid on StateNone
  Field                 * m_field_states[ MaximumFieldStates ];
  const DimensionTraits * m_dim_traits[ MaximumFieldDimension ];
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class Traits1 , class Traits2 ,
          class Traits3 , class Traits4 ,
          class Traits5 , class Traits6 ,
          class Traits7 >
class Field : public FieldBase {
public:

  typedef Scalar  data_type ;
  typedef Traits1 dimension_traits_1 ;
  typedef Traits2 dimension_traits_2 ;
  typedef Traits3 dimension_traits_3 ;
  typedef Traits4 dimension_traits_4 ;
  typedef Traits5 dimension_traits_5 ;
  typedef Traits6 dimension_traits_6 ;
  typedef Traits7 dimension_traits_7 ;

  typedef phdmesh::Dimension<Traits1,Traits2,Traits3,Traits4,
                             Traits5,Traits6,Traits7> Dimension ;

/*
  typedef phdmesh::Dimension<Traits1,Traits2,Traits3,Traits4,
                             Traits5,Traits6,Traits7,
                             EntityDimension> BlockDimension ;
*/

  Dimension dimension( unsigned entity_type , const Part & part ) const
    {
      const FieldBase::Dim & d = FieldBase::dimension( entity_type , part );
      return Dimension( d.stride );
    }

  const Field & operator[]( FieldState state ) const
    { return static_cast<Field &>( * FieldBase::m_field_states[state] ); }

private:

  ~Field();
  Field();
  Field( const Field & );
  Field & operator = ( const Field & );
};


} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

