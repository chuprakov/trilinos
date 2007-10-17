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

#include <strings.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <util/FixedArray.hpp>
#include <mesh/Field.hpp>
#include <mesh/Part.hpp>
#include <mesh/Schema.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const char * field_state_name( FieldState s )
{
  static const char * name_list[] = {
    "StateNew" ,
    "StateOld" ,
    "StateNM2" ,
    "StateNM3" ,
    "StateNM4" ,
    "StateNM5" ,
    "ERROR" };

  unsigned i = s ;
  if ( StateNM4 < i ) { i = 6 ; }
  return name_list[i] ;
}

//----------------------------------------------------------------------

FieldDimension::FieldDimension() : m_part( NULL )
{
  const unsigned zero = 0 ;
  Copy<MaxInfo>( m_info , zero );
}

FieldDimension::FieldDimension( const FieldDimension & rhs )
  : m_part( rhs.m_part )
{
  Copy<MaxInfo>( m_info , rhs.m_info );
}

FieldDimension & FieldDimension::operator = ( const FieldDimension & rhs )
{
  m_part = rhs.m_part ;
  Copy<MaxInfo>( m_info , rhs.m_info );
  return *this ;
}

FieldDimension::FieldDimension(
  const Part & p ,
  unsigned scalar_size ,
  unsigned n0 , unsigned n1 , unsigned n2 ,
  unsigned n3 , unsigned n4 , unsigned n5 ,
  unsigned n6 , unsigned n7 )
  : m_part( & p )
{
  enum { OK = StaticAssert< MaxDim == 8 >::OK };

  m_info[0] = 1 ;  m_info[1] = n0 ; m_info[2] = n1 ;
  m_info[3] = n2 ; m_info[4] = n3 ; m_info[5] = n4 ;
  m_info[6] = n5 ; m_info[7] = n6 ; m_info[8] = n7 ;
  m_info[7] = 0 ;  m_info[8] = 0 ;  m_info[9] = 0 ;

  for ( unsigned i = 0 ; i < MaxDim ; ++i ) {
    const unsigned n = i + 1 ;
    if ( m_info[n] *= m_info[i] ) {
      m_info[ I_Length ] = m_info[n] ;
      m_info[ I_Size ]   = m_info[n] * scalar_size ;
      m_info[ I_NDim ]   = n ;
    }
  }
}

unsigned FieldDimension::dimension( unsigned * const d ) const
{
  for ( unsigned i = 0 ; i < m_info[ I_NDim ] ; ++i ) {
    d[i] = m_info[i+1] / m_info[i] ;
  }
  return m_info[ I_NDim ];
}

bool FieldDimension::indices( unsigned off , unsigned * const ind ) const
{
  const bool result = off < length();

  if ( result ) {
    unsigned n = off ;
    unsigned i = number_of_dimensions();
    while ( i ) {
      --i ;
      ind[i] = n / m_info[i] ;
      n %= m_info[i] ;
    }
  }

  return result ;
}

bool FieldDimension::compare_dimension( const FieldDimension & rhs ) const
{
  return Compare<MaxInfo>::equal( m_info , rhs.m_info );
}

//----------------------------------------------------------------------

namespace {

struct LessFieldDimension {
  LessFieldDimension() {}

  bool operator()( const FieldDimension & lhs , const Part & rhs ) const 
  {
    const unsigned L = lhs.part()->schema_ordinal();
    const unsigned R = rhs.schema_ordinal();
    return L < R ;
  }
};

std::vector<FieldDimension>::iterator
find( std::vector<FieldDimension> & v , const Part & p )
{
  LessFieldDimension cmp ;
  std::vector<FieldDimension>::iterator
    i = std::lower_bound( v.begin() , v.end() , p , cmp );
  if ( i != v.end() && i->part() != & p ) { i = v.end(); }
  return i ;
}

std::vector<FieldDimension>::const_iterator
find( const std::vector<FieldDimension> & v , const Part & p )
{
  LessFieldDimension cmp ;
  std::vector<FieldDimension>::const_iterator
    i = std::lower_bound( v.begin() , v.end() , p , cmp );
  if ( i != v.end() && i->part() != & p ) { i = v.end(); }
  return i ;
}

void insert( std::vector<FieldDimension> & v , const FieldDimension & d )
{
  LessFieldDimension cmp ;
  const Part & p = * d.part();

  std::vector<FieldDimension>::iterator
    i = std::lower_bound( v.begin() , v.end() , p , cmp );

  if ( i == v.end() || i->part() != d.part() ) { v.insert( i , d ); }
}


}

//----------------------------------------------------------------------

Field<void,0>::~Field()
{ }

Field<void,0>::Field(
  Schema &            arg_schema ,
  EntityType          arg_type ,
  const std::string & arg_name ,
  unsigned scalar_type ,
  unsigned number_of_dimensions ,
  unsigned number_of_states ,
  FieldState state )
: m_cset() ,
  m_name( arg_name ),
  m_schema( arg_schema ),
  m_schema_ordinal(0),
  m_entity_type( arg_type ) ,
  m_scalar_type( scalar_type ),
  m_num_dim( number_of_dimensions ) ,
  m_num_states( number_of_states ),
  m_this_state( state ),
  m_dim_map()
{
  Field<void,0> * const pzero = NULL ;
  Copy<MaximumFieldStates>( m_field_states , pzero );
}

const std::vector<FieldDimension> & Field<void,0>::dimension() const
{ return m_field_states[0]->m_dim_map ; }

std::vector<FieldDimension> & Field<void,0>::dimension()
{ return m_field_states[0]->m_dim_map ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

struct FieldLessName {
  bool operator()( const Field<void,0> * lhs , const std::string & rhs ) const
    {
      const char * const l_name = lhs->name().c_str();
      const char * const r_name = rhs.c_str();
      return strcasecmp( l_name , r_name ) < 0 ;
    }
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

std::ostream &
get_field_header( std::ostream      & arg_msg ,
                  bool                arg_declare ,
                  EntityType          arg_entity_type ,
                  const std::string & arg_name ,
                  unsigned            arg_scalar_type ,
                  unsigned            arg_num_dim ,
                  unsigned            arg_num_states ,
                  const char        * arg_required_by )
{
  if ( arg_declare ) { arg_msg << "phdmesh::Schema::declare_field" ; }
  else               { arg_msg << "phdmesh::Schema::get_field" ; }
  arg_msg << "< " ;
  arg_msg << NumericEnum<>::name( arg_scalar_type );
  arg_msg << " , " ;
  arg_msg << arg_num_dim ;
  arg_msg << " >( " ;
  arg_msg << entity_type_name( arg_entity_type );
  arg_msg << " , " ;
  arg_msg << arg_name ;
  if ( arg_declare )     { arg_msg << " , " << arg_num_states ; }
  if ( arg_required_by ) { arg_msg << " , " << arg_required_by ; }
  arg_msg << " )" ;
  return arg_msg ;
}

}

Field<void,0> *
Schema::get_field(
  bool                arg_declare ,
  EntityType          arg_entity_type ,
  const std::string & arg_name ,
  unsigned            arg_scalar_type ,
  unsigned            arg_num_dim ,
  unsigned            arg_num_states ,
  const char * required_by ) const
{
  //----------------------------------------
  // ERROR CHECKING of input parameters:

  {
    const bool bad_entity_type = Other < arg_entity_type ;
    const bool bad_num_states  = arg_num_states < 1 ||
                                 MaximumFieldStates < arg_num_states ;
  
    if ( bad_entity_type || ( arg_declare && bad_num_states ) ) {
      std::ostringstream msg ;

      get_field_header( msg , arg_declare , arg_entity_type , arg_name ,
                              arg_scalar_type , arg_num_dim , arg_num_states ,
                              required_by );

      msg << " FAILED WITH BAD ARGUMENT" ;

      throw std::invalid_argument( msg.str() );
    }
  }

  // END ERROR CHECKING
  //----------------------------------------
  // ACTUAL WORK:

  const std::vector< Field<void,0> * > & field_set =
    m_fields[ arg_entity_type ];

  const std::vector< Field<void,0> * >::const_iterator e = field_set.end();
  const std::vector< Field<void,0> * >::const_iterator b = field_set.begin();

  const std::vector< Field<void,0> * >::const_iterator j =
    std::lower_bound( b, e, arg_name, FieldLessName() );

  Field<void,0> * const field =
     ( j != e && (*j)->name() == arg_name ) ? *j : NULL ;

  //----------------------------------------
  // ERROR CHECKING for compatibility:

  if ( field != NULL ) {

    const bool bad_type    = arg_scalar_type != field->numeric_type_ordinal();
    const bool bad_num_dim = arg_num_dim     != field->number_of_dimensions();
    const bool bad_num_states = arg_num_states != field->number_of_states() ;

    if ( bad_type || bad_num_dim || ( arg_declare && bad_num_states ) ) {
      std::ostringstream msg ;

      get_field_header( msg , arg_declare , arg_entity_type , arg_name ,
                              arg_scalar_type , arg_num_dim , arg_num_states ,
                              required_by );

      msg << " FAILED WITH INCOMPATIBLE " ;
      msg << "Field<" ;
      msg << NumericEnum<>::name( field->numeric_type_ordinal() );
      msg << "," ;
      msg << field->number_of_dimensions();
      msg << ">[" << arg_num_states << "]" ;

      throw std::invalid_argument( msg.str() );
    }
  }
  else if ( required_by ) {
    std::ostringstream msg ;

    get_field_header( msg , arg_declare , arg_entity_type , arg_name ,
                            arg_scalar_type , arg_num_dim , arg_num_states ,
                            required_by );
    msg << " FIND FAILURE" ;
    throw std::runtime_error( msg.str() );
  }

  // END ERROR CHECKING
  //----------------------------------------

  return field ;
}

//----------------------------------------------------------------------

Field<void,0> &
Schema::declare_field(
  EntityType          arg_entity_type,
  const std::string & arg_name ,
  unsigned            arg_scalar_type ,
  unsigned            arg_num_dim ,
  unsigned            arg_num_states )
{
  static const char method[] = "phdmesh::Schema::declare_field" ;

  assert_not_committed( method );

  Field<void,0> * field = get_field( true ,
                                     arg_entity_type , arg_name ,
                                     arg_scalar_type , arg_num_dim ,
                                     arg_num_states , NULL );

  if ( field == NULL ) {

    std::vector< Field<void,0> * > & field_set = m_fields[ arg_entity_type ];

    Field<void,0> * f[ MaximumFieldStates ] ;

    std::string field_names[ MaximumFieldStates ] ;

    switch( arg_num_states ) {
    case 6 : field_names[5] = arg_name ; field_names[5].append("_nm5");
    case 5 : field_names[4] = arg_name ; field_names[4].append("_nm4");
    case 4 : field_names[3] = arg_name ; field_names[3].append("_nm3");
    case 3 : field_names[2] = arg_name ; field_names[2].append("_nm2");
             field_names[1] = arg_name ; field_names[1].append("_nm1");
             break ;
    case 2 : field_names[1] = arg_name ; field_names[1].append("_old");
             break ;
    }

    field_names[0] = arg_name ;

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {

      const std::vector< Field<void,0> * >::iterator e = field_set.end();
      const std::vector< Field<void,0> * >::iterator b = field_set.begin();

      const std::vector< Field<void,0> * >::iterator
        j = std::lower_bound( b, e, field_names[i], FieldLessName() );

      if ( j != e && (*j)->name() == field_names[i] ) {
        Field<void,0> & tmp = **j ;
        std::ostringstream msg ;
        msg << method ;
        msg << " FAILED due to name collision with " ;
        msg << "Field< " ;
        msg << NumericEnum<>::name( tmp.numeric_type_ordinal() );
        msg << " , " ;
        msg << tmp.number_of_dimensions();
        msg << " > " ;
        msg << tmp.name();
        msg << "( " ;
        msg << entity_type_name( tmp.entity_type() );
        msg << " " ;
        msg << " States[" ;
        msg << tmp.number_of_states();
        msg << "] )" ;
        throw std::runtime_error( msg.str() );
      }

      f[i] = new Field<void,0>( *this , arg_entity_type,
                                field_names[i] ,
                                arg_scalar_type, arg_num_dim,
                                arg_num_states , (FieldState) i );

      field_set.insert( j , f[i] );
    }

    field = f[0] ;

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      Field<void,0> & tmp = * f[i] ;
      for ( unsigned k = 0 ; k < arg_num_states ; ++k ) {
        tmp.m_field_states[k] = f[k] ;
      }
    }

    // Update ordinals to reflect the newly inserted fields
    {
      const std::vector< Field<void,0> *>::iterator e = field_set.end();
            std::vector< Field<void,0> *>::iterator j = field_set.begin();
      for ( unsigned k = 0 ; j != e ; ++j , ++k ) {
        Field<void,0> & tmp = **j ;
        tmp.m_schema_ordinal = k ;
        if ( tmp.m_field_states[0]->m_schema_ordinal + tmp.m_this_state != k ) {
          std::string msg ;
          msg.append( method );
          msg.append( " FAILED ORDERING OF FIELDS" );
          throw std::logic_error(msg);
        }
      }
    }
  }

  return *field ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void assert_field_dimension_compatible(
  const char * const method ,
  const char * const field_name ,
  const FieldDimension & dim_a ,
  const FieldDimension & dim_b )
{
  if ( ! dim_a.compare_dimension( dim_b ) ) {
    std::ostringstream msg ;
    msg << method ;
    msg << "[ " ;
    msg << field_name ;
    msg << " ]" ;
    msg << " FAILED: incompatible dimensions : { " ;

    dim_a.part()->name();
    for ( unsigned i = 0 ; i < MaximumFieldDimension ; ++i ) {
      msg << " " << dim_a[i] ;
    }
    msg << " } != { " ;
    dim_b.part()->name();
    for ( unsigned i = 0 ; i < MaximumFieldDimension ; ++i ) {
      msg << " " << dim_b[i] ;
    }
    msg << " }" ;

    throw std::runtime_error( msg.str() );
  }
}

}

//----------------------------------------------------------------------
// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void Field<void,0>::set_dimension(
  const Part & p , unsigned n0 , unsigned n1 , unsigned n2 ,
                   unsigned n3 , unsigned n4 , unsigned n5 ,
                   unsigned n6 , unsigned n7 )
{
  static const char method[] = "phdmesh::Field::set_dimension" ;

  m_schema.assert_not_committed( method );

  FieldDimension dim( p , NumericEnum<>::size( m_scalar_type ) ,
                      n0 , n1 , n2 , n3 , n4 , n5 , n6 , n7 );

  std::vector<FieldDimension> & dim_map = dimension();

  //----------------------------------------
  // ERROR CHECKING:
  {
    const unsigned type_size  = NumericEnum<>::size( m_scalar_type );
    const unsigned dim_length = dim.length();
    const unsigned dim_size   = dim.size();
    const unsigned dim_number = dim.number_of_dimensions();
    const bool bad_num_dim    = m_num_dim != dim_number ;
    const bool bad_type_size  = dim_length * type_size != dim_size ;

    if ( bad_num_dim || bad_type_size ) {
      std::ostringstream msg ;
      msg << method ;
      msg << " FAILED WITH INCOMPATIBLE DIMENSION : " ;
      msg << entity_type_name( m_entity_type );
      msg << " Field< " ;
      msg << NumericEnum<>::name( m_scalar_type );
      msg << "(" << type_size ;
      msg << ") , " << m_num_dim ;
      msg << " > " ;
      msg << m_name ;
      if ( bad_num_dim ) {
        msg << " #Dimensions = " ;
        msg << dim_number ;
      }
      else {
        msg << " TypeSize(" ;
        msg << ( dim_length ? ( dim_size / dim_length ) : 0 );
        msg << ")" ;
      }
      throw std::runtime_error( msg.str() );
    }
  }
  // END ERROR CHECKING
  //----------------------------------------
  // Retrieve dimension for this part, also retrieves for supersets.

  const FieldDimension & old_dim = dimension( p );

  if ( old_dim.length() ) { // This part or a superset already exists
    assert_field_dimension_compatible( method , m_name.c_str(),
                                       dim , old_dim );
  }
  else {

    PartSet::const_iterator ip ;

    // All subsets must be compatible:

    for ( ip = p.subsets().begin() ; ip != p.subsets().end() ; ++ip ) {
      const FieldDimension & tmp_dim = dimension( **ip );
      if ( tmp_dim.length() ) {
        assert_field_dimension_compatible( method, m_name.c_str(),
                                           dim, tmp_dim );
      }
    }

    // If subset already present then remove it:

    for ( ip = p.subsets().begin() ; ip != p.subsets().end() ; ++ip ) {
      std::vector<FieldDimension>::iterator itmp = find( dim_map , **ip );
      if ( itmp != dim_map.end() ) { dim_map.erase( itmp ); }
    }

    insert( dim_map , dim );
  }
}

//----------------------------------------------------------------------
// If a part and one of its subset parts show up in the dimension map
// verify compatibility of dimensions and delete the subset part.

void Field<void,0>::clean_dimension()
{
  typedef std::vector<FieldDimension>::difference_type d_type ;

  static const char method[] = "phdmesh::Field::clearn_dimension" ;

  m_schema.assert_not_committed( method );

  std::vector<FieldDimension> & dim_map = dimension();

  std::vector<FieldDimension>::iterator i = dim_map.begin();

  for ( ; i != dim_map.end() ; ++i ) {
    FieldDimension & dim = *i ;
    const PartSet  & sub = dim.part()->subsets();

    std::vector<FieldDimension>::iterator j = dim_map.begin();

    while ( j != dim_map.end() ) {
      std::vector<FieldDimension>::iterator k = j ; ++j ;
      if ( i != k ) {
        if ( contain( sub , * k->part() ) ) {
          assert_field_dimension_compatible( method, m_name.c_str(),
                                             dim, *k );
          dim_map.erase( k );
        }
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// This part or any superset of this part

const FieldDimension &
Field<void,0>::dimension( const Part & part ) const
{
  static const FieldDimension empty ;

  const std::vector<FieldDimension> & dim_map = dimension();

  const std::vector<FieldDimension>::const_iterator ie = dim_map.end() ;

  std::vector<FieldDimension>::const_iterator i = find( dim_map , part );

  const PartSet::const_iterator ipe = part.supersets().end();
        PartSet::const_iterator ip  = part.supersets().begin() ;

  for ( ; i == ie && ip != ipe ; ++ip ) {
    const Part & p = **ip ;
    i = find( dim_map , p );
  }

  return i == ie ? empty : *i ;
}

unsigned Field<void,0>::max_length() const
{
  unsigned max = 0 ;

  const std::vector<FieldDimension> & dim_map = dimension();

  const std::vector<FieldDimension>::const_iterator ie = dim_map.end();
        std::vector<FieldDimension>::const_iterator i  = dim_map.begin();

  for ( ; i != ie ; ++i ) {
    const unsigned len = i->length();
    if ( max < len ) { max = len ; }
  }

  return max ;
}

unsigned Field<void,0>::max_size() const
{
  unsigned max = 0 ;

  const std::vector<FieldDimension> & dim_map = dimension();

  const std::vector<FieldDimension>::const_iterator ie = dim_map.end();
        std::vector<FieldDimension>::const_iterator i  = dim_map.begin();

  for ( ; i != ie ; ++i ) {
    const unsigned siz = i->size();
    if ( max < siz ) { max = siz ; }
  }

  return max ;
}

//----------------------------------------------------------------------

void Field<void,0>::assert_validity( unsigned ftype ,
                                     unsigned ndim ,
                                     unsigned state ) const
{
  enum { vtype = NumericEnum<void>::value };

  const bool bad_type  = vtype != ftype && ftype != m_scalar_type ;
  const bool bad_ndim  = 0     != ndim  && ndim  != m_num_dim ;
  const bool bad_state = m_num_states <= state ;

  if ( bad_type || bad_ndim || bad_state ) {
    std::ostringstream msg ;
    msg << "phdmesh::" ;
    msg << entity_type_name( m_entity_type );
    msg << " Field< " ;
    msg << NumericEnum<>::name( m_scalar_type );
    msg << " , " ;
    msg << m_num_dim ;
    msg << " > " ;
    msg << m_name ;
    if ( bad_state ) {
      msg << " [ #States = " ;
      msg << m_num_states ;
      msg << " ] " ;
    }
    msg << " ERROR INVALID CAST TO Field<" ;
    msg << NumericEnum<>::name( ftype );
    msg << " , " ;
    msg << ndim ;
    msg << " >" ;
    if ( bad_state ) {
      msg << " [ State = " ;
      msg << field_state_name( FieldState(state) );
      msg << "(" << state << ")" ;
    }
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------

CSet & Field<void,0>::cset_update()
{
  static const char method[] = "phdmesh::Schema::cset_update" ;

  m_schema.assert_not_committed( method );

  return m_cset ;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const Field<void,0> & f )
{
  s << "Field< " ;
  s << NumericEnum<>::name( f.numeric_type_ordinal() );
  s << " , " ;
  s << f.number_of_dimensions();
  s << " >( " ;
  s << entity_type_name( f.entity_type() );
  s << " , " ;
  s << f.name() ;
  if ( 1 < f.number_of_states() ) {
    s << "[ " ;
    s << field_state_name( f.state() );
    s << " ]" ;
  }
  s << " )" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b , const Field<void,0> & f )
{
  const std::vector<FieldDimension> & dim_map = f.dimension();
  s << f ;
  s << " {" ;
  for ( std::vector<FieldDimension>::const_iterator
        i = dim_map.begin() ; i != dim_map.end() ; ++i ) {
    const FieldDimension & d = *i ;
    s << std::endl << b << "  " ;
    s << d.part()->name();
    s << " ( " << d[0] ;
    const unsigned n = d.number_of_dimensions();
    for ( unsigned j = 1 ; j < n ; ++j ) {
      s << " , " << d[j] ;
    }
    s << " )" ;
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

