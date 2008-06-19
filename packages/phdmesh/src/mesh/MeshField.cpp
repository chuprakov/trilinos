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
#include <mesh/Types.hpp>
#include <mesh/Field.hpp>
#include <mesh/Part.hpp>
#include <mesh/Schema.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

const DimensionTraits * EntityDimension::descriptor()
{ static const EntityDimension self ; return & self ; }

const char * EntityDimension::name() const
{ static const char n[] = "EntityDimension" ; return n ; }

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

namespace {

struct LessFieldDimension {
  LessFieldDimension() {}

  bool operator()( const FieldBase::Dim & lhs , const FieldBase::Dim & rhs )
  {
    unsigned L = lhs.type ;
    unsigned R = rhs.type ;
    if ( L == R ) {
      L = lhs.part->schema_ordinal();
      R = rhs.part->schema_ordinal();
    }
    return L < R ;
  }
};

std::vector<FieldBase::Dim>::const_iterator
find( const std::vector<FieldBase::Dim> & v , const FieldBase::Dim & p )
{
  LessFieldDimension cmp ;
  std::vector<FieldBase::Dim>::const_iterator
    i = std::lower_bound( v.begin() , v.end() , p , cmp );
  if ( i != v.end() && ( i->part != p.part || i->type != p.type ) ) {
    i = v.end();
  }
  return i ;
}

void insert( std::vector<FieldBase::Dim> & v , const FieldBase::Dim & d )
{
  LessFieldDimension cmp ;

  std::vector<FieldBase::Dim>::iterator
    i = std::lower_bound( v.begin() , v.end() , d , cmp );

  if ( i == v.end() || i->part != d.part || i->type != d.type ) {
    v.insert( i , d );
  }
}


}

//----------------------------------------------------------------------

FieldBase::~Field()
{ }

FieldBase::Field(
  Schema &            arg_schema ,
  const std::string & arg_name ,
  unsigned scalar_type ,
  const DimensionTraits * t1 ,
  const DimensionTraits * t2 ,
  const DimensionTraits * t3 ,
  const DimensionTraits * t4 ,
  const DimensionTraits * t5 ,
  const DimensionTraits * t6 ,
  const DimensionTraits * t7 ,
  unsigned number_of_states ,
  FieldState state )
: m_cset() ,
  m_name( arg_name ),
  m_schema( arg_schema ),
  m_schema_ordinal(0),
  m_scalar_type( scalar_type ),
  m_num_dim( ( ! t1 ? 0 :
             ( ! t2 ? 1 :
             ( ! t3 ? 2 :
             ( ! t4 ? 3 :
             ( ! t5 ? 4 :
             ( ! t6 ? 5 :
             ( ! t7 ? 6 : 7 )))))))),
  m_num_states( number_of_states ),
  m_this_state( state ),
  m_dim_map()
{
  FieldBase * const pzero = NULL ;
  Copy<MaximumFieldStates>( m_field_states , pzero );
  m_dim_traits[0] = t1 ;
  m_dim_traits[1] = t2 ;
  m_dim_traits[2] = t3 ;
  m_dim_traits[3] = t4 ;
  m_dim_traits[4] = t5 ;
  m_dim_traits[5] = t6 ;
  m_dim_traits[6] = t7 ;
}

const std::vector<FieldBase::Dim> & FieldBase::dimension() const
{ return m_field_states[0]->m_dim_map ; }

std::vector<FieldBase::Dim> & FieldBase::dimension()
{ return m_field_states[0]->m_dim_map ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

struct FieldLessName {
  bool operator()( const FieldBase * lhs , const std::string & rhs ) const
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

void
print_field_type( std::ostream          & arg_msg ,
                  unsigned                arg_scalar_type ,
                  const DimensionTraits * arg_t1 ,
                  const DimensionTraits * arg_t2 ,
                  const DimensionTraits * arg_t3 ,
                  const DimensionTraits * arg_t4 ,
                  const DimensionTraits * arg_t5 ,
                  const DimensionTraits * arg_t6 ,
                  const DimensionTraits * arg_t7 )
{
  arg_msg << "Field< " ;
  arg_msg << NumericEnum<>::name( arg_scalar_type );
  if ( arg_t1 ) { arg_msg << " , " << arg_t1->name();
  if ( arg_t2 ) { arg_msg << " , " << arg_t2->name();
  if ( arg_t3 ) { arg_msg << " , " << arg_t3->name();
  if ( arg_t4 ) { arg_msg << " , " << arg_t4->name();
  if ( arg_t5 ) { arg_msg << " , " << arg_t5->name();
  if ( arg_t6 ) { arg_msg << " , " << arg_t6->name();
  if ( arg_t7 ) { arg_msg << " , " << arg_t7->name();
  } } } } } } }
  arg_msg << " >" ;
}

}

FieldBase *
Schema::get_field_base(
  const std::string & arg_name ,
  unsigned          arg_scalar_type ,
  const DimensionTraits * arg_t1 ,
  const DimensionTraits * arg_t2 ,
  const DimensionTraits * arg_t3 ,
  const DimensionTraits * arg_t4 ,
  const DimensionTraits * arg_t5 ,
  const DimensionTraits * arg_t6 ,
  const DimensionTraits * arg_t7 ,
  int arg_number_states ,
  const char * arg_required_by ) const
{
  static const char declare_method[] = "phdmesh::Schema::declare_field" ;
  static const char get_method[]     = "phdmesh::Schema::get_field" ;

  // Potential error conditions:

  bool not_found         = false ;
  bool bad_scalar_type   = false ;
  bool bad_dimension     = false ;
  bool bad_number_states = MaximumFieldStates < arg_number_states ;

  // Find the field by name:

  const std::vector< FieldBase * >::const_iterator e = m_fields.end();
  const std::vector< FieldBase * >::const_iterator b = m_fields.begin();

  const std::vector< FieldBase * >::const_iterator j =
    std::lower_bound( b, e, arg_name, FieldLessName() );

  FieldBase * const field =
     ( j != e && (*j)->name() == arg_name ) ? *j : NULL ;

  // If found check for compatibility:

  if ( field != NULL ) {

    bad_scalar_type = arg_scalar_type != field->numeric_type_ordinal();

    bad_dimension = arg_t1 != field->m_dim_traits[0] ||
                    arg_t2 != field->m_dim_traits[1] ||
                    arg_t3 != field->m_dim_traits[2] ||
                    arg_t4 != field->m_dim_traits[3] ||
                    arg_t5 != field->m_dim_traits[4] ||
                    arg_t6 != field->m_dim_traits[5] ||
                    arg_t7 != field->m_dim_traits[6] ;

    bad_number_states =
      ( 0 <= arg_number_states ) &&
      ( arg_number_states != (int) field->m_num_states );
  }
  else if ( arg_required_by ) {
    not_found = true ;
  }

  if ( bad_scalar_type || bad_dimension || bad_number_states || not_found ) {

    const char * const method =
      arg_number_states < 0 ? get_method : declare_method ;

    std::ostringstream msg ;

    msg << method << "< " ;
    print_field_type( msg , arg_scalar_type ,
                            arg_t1 , arg_t2 , arg_t3 ,
                            arg_t4 , arg_t5 , arg_t6 , arg_t7 );
    msg << " >( \"" << arg_name << "\"" ;

    if ( 0 <= arg_number_states ) {
      msg << " , " << arg_number_states ;
      if ( bad_number_states ) { msg << " <= NUMBER_OF_STATES IS BAD" ; }
    }

    if ( arg_required_by ) {
      msg << " , " << arg_required_by << " <= REQUIRED_BY" ;
      if ( not_found ) { msg << " BUT NOT FOUND" ; }
    }

    msg << " )" ;
 
    if ( field != NULL ) {
      msg << " FOUND INCOMPATIBLE " ;
      print_field_type( msg , field->numeric_type_ordinal() ,
                              field->m_dim_traits[0] ,
                              field->m_dim_traits[1] ,
                              field->m_dim_traits[2] ,
                              field->m_dim_traits[3] ,
                              field->m_dim_traits[4] ,
                              field->m_dim_traits[5] ,
                              field->m_dim_traits[6] );
      msg << "( " << field->m_num_states << " )" ;
    }
    throw std::runtime_error( msg.str() );
  }

  return field ;
}

//----------------------------------------------------------------------

FieldBase &
Schema::declare_field_base(
  const std::string & arg_name ,
  unsigned            arg_scalar_type ,
  const DimensionTraits * arg_t1 ,
  const DimensionTraits * arg_t2 ,
  const DimensionTraits * arg_t3 ,
  const DimensionTraits * arg_t4 ,
  const DimensionTraits * arg_t5 ,
  const DimensionTraits * arg_t6 ,
  const DimensionTraits * arg_t7 ,
  unsigned            arg_num_states )
{
  static const char method[] = "phdmesh::Schema::declare_field" ;

  static const char reserved_state_suffix[6][5] = {
    "_OLD" , "_NM1" , "_NM2" , "_NM3" , "_NM4" , "_NM5" };

  assert_not_committed( method );

  // Check that the name does not have a reserved suffix

  for ( unsigned i = 0 ; i < 6 ; ++i ) {
    const int len = arg_name.size() - 4 ;
    if ( 0 <= len &&
         ! strcasecmp( arg_name.c_str() + len , reserved_state_suffix[i] ) ) {
      std::ostringstream msg ;
      msg << method << "< " ;
      print_field_type( msg , arg_scalar_type ,
                              arg_t1 , arg_t2 , arg_t3 ,
                              arg_t4 , arg_t5 , arg_t6 , arg_t7 );
      msg << " >( \"" << arg_name ;
      msg << "\" <= HAS RESERVED STATE SUFFIX \"" ;
      msg << reserved_state_suffix[i] ;
      msg << " )" ;
      throw std::runtime_error( msg.str() );
    }
  }

  // Get field and if found verify compatibility

  FieldBase * field = get_field_base( arg_name ,
                                      arg_scalar_type ,
                                      arg_t1 , arg_t2 , arg_t3 ,
                                      arg_t4 , arg_t5 , arg_t6 , arg_t7 ,
                                      arg_num_states , NULL );

  if ( field == NULL ) {

    FieldBase * f[ MaximumFieldStates ] ;

    std::string field_names[ MaximumFieldStates ] ;

    field_names[0] = arg_name ;

    if ( 2 == arg_num_states ) {
      field_names[1] = arg_name ;
      field_names[1].append( reserved_state_suffix[0] );
    }
    else {
      for ( unsigned i = 1 ; i < arg_num_states ; ++i ) {
        field_names[i] = arg_name ;
        field_names[i].append( reserved_state_suffix[i] );
      }
    }

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {

      const std::vector< FieldBase * >::iterator e = m_fields.end();
      const std::vector< FieldBase * >::iterator b = m_fields.begin();

      const std::vector< FieldBase * >::iterator
        j = std::lower_bound( b, e, field_names[i], FieldLessName() );

      if ( j != e && (*j)->name() == field_names[i] ) {
        std::string msg( method );
        msg.append(" CATASTROPHIC INTERNAL LOGIC ERROR FOR NAMES" );
        throw std::logic_error( msg );
      }

      f[i] = new FieldBase( *this ,
                            field_names[i] ,
                            arg_scalar_type,
                            arg_t1 , arg_t2 , arg_t3 ,
                            arg_t4 , arg_t5 , arg_t6 , arg_t7 ,
                            arg_num_states , (FieldState) i );

      m_fields.insert( j , f[i] );
    }

    field = f[0] ;

    for ( unsigned i = 0 ; i < arg_num_states ; ++i ) {
      FieldBase & tmp = * f[i] ;
      for ( unsigned k = 0 ; k < arg_num_states ; ++k ) {
        tmp.m_field_states[k] = f[k] ;
      }
    }

    // Update ordinals to reflect the newly inserted fields
    {
      const std::vector< FieldBase *>::iterator e = m_fields.end();
            std::vector< FieldBase *>::iterator j = m_fields.begin();
      for ( unsigned k = 0 ; j != e ; ++j , ++k ) {
        FieldBase & tmp = **j ;
        tmp.m_schema_ordinal = k ;
        if ( tmp.m_field_states[0]->m_schema_ordinal + tmp.m_this_state != k ) {
          std::string msg( method );
          msg.append(" CATASTROPHIC INTERNAL LOGIC ERROR FOR ORDINALS" );
          throw std::logic_error(msg);
        }
      }
    }
  }

  return *field ;
}

//----------------------------------------------------------------------

void Schema::declare_field_relation(
  FieldBase & pointer_field ,
  relation_stencil_ptr stencil ,
  FieldBase & referenced_field )
{
  static const char method[] = "phdmesh::Schema::declare_field_relation" ;

  static const int offset = NumericEnum<void*>::value -
                            NumericEnum<void>::value ;

  if ( referenced_field.numeric_type_ordinal() + offset !=
       pointer_field.numeric_type_ordinal() ) {

    std::string msg( method );
    msg.append( ": FAILED, " );
    msg.append( pointer_field.name() );
    msg.append( " and " );
    msg.append( referenced_field.name() );
    msg.append( " have incompatible numerical types." );
    throw std::invalid_argument( msg );
  }

  FieldRelation tmp ;
  tmp.m_root   = & pointer_field ;
  tmp.m_target = & referenced_field ;
  tmp.m_function = stencil ;

  m_field_relations.push_back( tmp );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void print_field_dim( std::ostream & msg ,
                      const unsigned n ,
                      const FieldBase::Dim & d )
{
  msg << "( " << entity_type_name( d.type );
  msg << " , Part[" ;
  msg << d.part->name();
  msg << "] , " << d.stride[0] ;
  for ( unsigned j = 1 ; j < n ; ++j ) {
    if ( d.stride[j-1] ) {
      msg << " , " << ( d.stride[j] / d.stride[j-1] );
    }
    else {
      msg << " , ERROR " ;
    }
  }
  msg << ")" ;
}

void assert_field_dimension_compatible(
  const char * const method ,
  const FieldBase & field ,
  const FieldBase::Dim & a ,
  const FieldBase::Dim & b )
{
  if ( Compare< MaximumFieldDimension >::not_equal( a.stride , b.stride ) ) {
    std::ostringstream msg ;
    msg << method << " FOUND INCOMPATIBLE SIZES FOR " ;
    print_field_type( msg , field.numeric_type_ordinal() ,
                            field.dimension_traits()[0] ,
                            field.dimension_traits()[1] ,
                            field.dimension_traits()[2] ,
                            field.dimension_traits()[3] ,
                            field.dimension_traits()[4] ,
                            field.dimension_traits()[5] ,
                            field.dimension_traits()[6] );
    msg << "[" << field.name() << "] " ;

    print_field_dim( msg , field.number_of_dimensions() , a );
    msg << " INCOMPATIBLE WITH " ;
    print_field_dim( msg , field.number_of_dimensions() , b );
    throw std::runtime_error( msg.str() );
  }
}

}

//----------------------------------------------------------------------
// Setting the dimension for one field sets the dimension
// for the corresponding fields of the FieldState array.
// If subset exists then replace it.
// If exists or superset exists then do nothing.

void Schema::declare_field_stride(
  FieldBase      & arg_field ,
  EntityType       arg_entity_type ,
  const Part     & arg_part ,
  const unsigned * arg_stride )
{
  static const char method[] = "phdmesh::Schema::declare_field_dimension" ;

  assert_not_committed( method );
  assert_same_schema( method , arg_field.m_schema );
  assert_same_schema( method , arg_part.m_schema );

  FieldBase::Dim tmp( arg_entity_type , arg_part );
  {
    unsigned j ;
    for ( j = 0 ; j < arg_field.m_num_dim ; ++j ) {
      tmp.stride[j] = arg_stride[j] ;
    }
    for ( ; j < MaximumFieldDimension ; ++j ) {
      tmp.stride[j] = 0 ;
    }
  }

  // If a dimension is defined for a superset of the part
  // then that dimension must be compatible and this
  // declaration is redundant.

  bool redundant  = false ;

  std::vector<FieldBase::Dim> & dim_map =
    arg_field.m_field_states[0]->m_dim_map ;

  std::vector<FieldBase::Dim>::iterator i ;

  for ( i = dim_map.begin() ; i != dim_map.end() ; ) {
    bool remove = false ;

    if ( arg_entity_type == i->type ) {

      const Part & part = * i->part ;

      if ( ( arg_part == part ) || ( contain( arg_part.subsets() , part ) ) ) {
        redundant = true ;
      }
      else if ( contain( arg_part.supersets() , part ) ) {
        remove = true ;
      }

      if ( redundant || remove ) {
        assert_field_dimension_compatible( method , arg_field , tmp , *i );
      }
    }

    if ( remove ) {
      i = dim_map.erase( i );
    }
    else {
      ++i ;
    }
  }

  if ( ! redundant ) { insert( dim_map , tmp ); }
}

//----------------------------------------------------------------------
// If a part and one of its subset parts show up in the dimension map
// verify compatibility of dimensions and delete the subset part.

void Schema::clean_field_dimension()
{
  static const char method[] = "phdmesh::Schema::clean_field_dimension" ;
  const int zero = 0 ;

  for ( std::vector<FieldBase *>::iterator
        f = m_fields.begin() ; f != m_fields.end() ; ++f ) {
    FieldBase & field = **f ;

    std::vector<FieldBase::Dim> & dim_map = field.m_dim_map ;
    std::vector<int> flag( dim_map.size() , zero );

    for ( size_t i = 0 ; i < dim_map.size() ; ++i ) {
      const FieldBase::Dim & dim = dim_map[i] ;
      const PartSet  & sub = dim_map[i].part->subsets();

      for ( size_t j = 0 ; j < dim_map.size() ; ++j ) {
        if ( i != j &&
             dim.type == dim_map[j].type &&
             contain( sub , * dim_map[j].part ) ) {
          assert_field_dimension_compatible( method, field, dim, dim_map[j] );
          flag[j] = 1 ;
        }
      }
    }

    for ( size_t i = dim_map.size() ; i-- ; ) {
      if ( flag[i] ) {
        std::vector<FieldBase::Dim>::iterator j = dim_map.begin();
        std::advance( j , i );
        dim_map.erase( j );
      }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// This part or any superset of this part

const FieldBase::Dim &
FieldBase::dimension( EntityType etype , const Part & part ) const
{
  static const FieldBase::Dim empty ;

  FieldBase::Dim tmp( etype , part );

  const std::vector<FieldBase::Dim> & dim_map = dimension();

  const std::vector<FieldBase::Dim>::const_iterator ie = dim_map.end() ;

  std::vector<FieldBase::Dim>::const_iterator i = find( dim_map , tmp );

  const PartSet::const_iterator ipe = part.supersets().end();
        PartSet::const_iterator ip  = part.supersets().begin() ;

  for ( ; i == ie && ip != ipe ; ++ip ) {
    tmp.part = *ip ;
    i = find( dim_map , tmp );
  }

  return i == ie ? empty : *i ;
}

unsigned FieldBase::max_size( EntityType entity_type ) const
{
  unsigned max = 0 ;

  const std::vector<FieldBase::Dim> & dim_map = dimension();

  const std::vector<FieldBase::Dim>::const_iterator ie = dim_map.end();
        std::vector<FieldBase::Dim>::const_iterator i  = dim_map.begin();

  for ( ; i != ie ; ++i ) {
    if ( i->type == entity_type ) {
      const unsigned len = m_num_dim ? i->stride[ m_num_dim - 1 ] : 1 ;
      if ( max < len ) { max = len ; }
    }
  }

  return max ;
}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const FieldBase & field )
{
  print_field_type( s , field.numeric_type_ordinal() ,
                        field.dimension_traits()[0] ,
                        field.dimension_traits()[1] ,
                        field.dimension_traits()[2] ,
                        field.dimension_traits()[3] ,
                        field.dimension_traits()[4] ,
                        field.dimension_traits()[5] ,
                        field.dimension_traits()[6] );
  s << "[ \"" ;
  s << field.name() ;
  s << " \" , " ;
  s << field.number_of_states();
  s << " ]" ;
  return s ;
}

std::ostream & print( std::ostream & s ,
                      const char * const b ,
                      const FieldBase & field )
{
  const std::vector<FieldBase::Dim> & dim_map = field.dimension();
  s << field ;
  s << " {" ;
  for ( std::vector<FieldBase::Dim>::const_iterator
        i = dim_map.begin() ; i != dim_map.end() ; ++i ) {
    s << std::endl << b << "  " ;
    print_field_dim( s , field.number_of_dimensions() , *i );
  }
  s << std::endl << b << "}" ;
  return s ;
}

//----------------------------------------------------------------------

} // namespace phdmesh

