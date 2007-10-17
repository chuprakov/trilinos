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
#include <util/CSet.hpp>

#include <mesh/Types.hpp>
#include <mesh/FieldDim.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

/** if ( #states == 1 ) "Field<T,N>( entity , name )"
 *  else                "Field<T,N>( entity , name[ state ] )"
 */
std::ostream & operator << ( std::ostream & , const Field<void,0> & );

/** Print field and field dimension map entries on new lines. */
std::ostream & print( std::ostream & ,
                      const char * const , const Field<void,0> & );

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
class Field<void,0> {
public:

  Schema & schema() const { return m_schema ; }

  /** Volatile until schema is committed. */
  unsigned schema_ordinal() const { return m_schema_ordinal ; }

  EntityType entity_type() const { return m_entity_type ; }

  const std::string & name() const { return m_name ; }

  template<typename NumType> bool type_is() const
    { return m_scalar_type == NumericEnum<NumType>::value ; }

  unsigned numeric_type_ordinal() const { return m_scalar_type ; }

  unsigned number_of_dimensions() const { return m_num_dim ; }

  unsigned number_of_states() const { return m_num_states ; }

  FieldState state() const { return m_this_state ; }

  unsigned max_length() const ;

  unsigned max_size() const ;

  /** Volatile until the schema is committed */
  const FieldDimension & dimension( const Part & ) const ;

  /** Volatile until the schema is committed */
  const std::vector<FieldDimension> & dimension() const ;

  /** The schema cannot be committed.
   *  The number of non-zero arguments must match the number of dimensions.
   */
  void set_dimension( const Part & , unsigned n0 ,
                                     unsigned n1 = 0 , unsigned n2 = 0 ,
                                     unsigned n3 = 0 , unsigned n4 = 0 ,
                                     unsigned n5 = 0 , unsigned n6 = 0 ,
                                     unsigned n7 = 0 );

  /** Verify and clean part->dimension mapping */
  void clean_dimension();

  //----------------------------------------

  const CSet & cset_query() const { return m_cset ; }
        CSet & cset_update();

  //----------------------------------------
  /** Checked conversion to field of the given specification. */
  template<typename T, unsigned NDim>
    const Field<T,NDim> & field( FieldState s = StateNew ) const
      {
        assert_validity( NumericEnum<T>::value, NDim, (unsigned) s );
        return static_cast<const Field<T,NDim> &>( * m_field_states[s]  );
      }

private:

  friend class Schema ;

  void assert_validity( unsigned , unsigned , unsigned ) const ;

  ~Field();

  Field();
  Field( const Field<void,0> & );
  Field<void,0> & operator = ( const Field<void,0> & );

  Field( Schema & ,
         EntityType ,
         const std::string & ,
         unsigned scalar_type ,
         unsigned number_of_dimensions ,
         unsigned number_of_states ,
         FieldState );

  std::vector<FieldDimension> & dimension();

  CSet        m_cset ;
  std::string m_name ;
  Schema    & m_schema ;         // Mesh schema in which this field resides
  unsigned    m_schema_ordinal ; // Ordinal in the field set
  EntityType  m_entity_type ;    // Type of mesh entities
  unsigned    m_scalar_type ;    // Ordinal in FieldTypes
  unsigned    m_num_dim ;        // Number of dimensions
  unsigned    m_num_states ;     // Number of states of this field
  FieldState  m_this_state ;     // Field state of this field
  std::vector<FieldDimension> m_dim_map ; // Only valid on StateNone
  Field<void,0> * m_field_states[ MaximumFieldStates ];
};

/** Specify fundamental type and number of dimensions of a field */

template<typename T, unsigned NDim>
class Field : public Field<void,0> {
private:
  enum { OkType = NumericEnum<T>::OK };
  enum { OkDim  = StaticAssert< NDim <= MaximumFieldDimension >::OK };
public:

  const Field<T,NDim> & operator[]( FieldState state ) const
    { return Field<void,0>::field<T,NDim>( state ); }

private:
  ~Field();
  Field();
  Field( const Field<T,NDim> & );
  Field<T,NDim> operator = ( const Field<T,NDim> & );
};


} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

