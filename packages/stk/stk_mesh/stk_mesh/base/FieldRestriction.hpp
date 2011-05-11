/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_baseImpl_FieldRestriction_hpp
#define stk_mesh_baseImpl_FieldRestriction_hpp

#include <vector>
#include <Shards_Array.hpp>
#include <stk_mesh/base/EntityKey.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/util/SimpleArrayOps.hpp>

namespace stk {
namespace mesh {

/**
 * A field restrictions is one of the three fundamental components of a field specification
 * (see Field.hpp for a full discusssion); it defines a set of entities that have a field.
 *
 * This class encapsulates a minimal set of data for a field restriction. The API for
 * declaring field restrictions is in MetaData.hpp.
 */
class FieldRestriction {
  public:

  typedef shards::array_traits::int_t size_type ;

  FieldRestriction() : m_key() {
    Copy<MaximumFieldDimension>( m_stride , size_type(0) );
  }

  FieldRestriction( const FieldRestriction & rhs ): m_key( rhs.m_key ) {
    Copy< MaximumFieldDimension >( m_stride , rhs.m_stride );
  }

  FieldRestriction & operator = ( const FieldRestriction & rhs ) {
    m_key = rhs.m_key ;
    Copy< MaximumFieldDimension >( m_stride , rhs.m_stride );
    return *this ;
  }

  FieldRestriction( EntityRank input_rank , PartOrdinal input_ordinal)
    : m_key( input_rank, input_ordinal )
  {
    Copy< MaximumFieldDimension >( m_stride , size_type(0) );
  }

  EntityRank rank()    const { return entity_rank( m_key ); }
  PartOrdinal ordinal() const { return entity_id( m_key ); }

  size_type & stride_nonconst( Ordinal index ) { return m_stride[index]; }
  const size_type & stride( Ordinal index ) const { return m_stride[index]; }

  size_type dimension() const { return m_stride[0]; }

  bool operator < ( const FieldRestriction & rhs ) const { return this->m_key < rhs.m_key; }
  bool operator == ( const FieldRestriction & rhs ) const { return this->m_key == rhs.m_key; }
  bool operator != ( const FieldRestriction & rhs ) const { return this->m_key != rhs.m_key; }

  bool not_equal_stride( const FieldRestriction & rhs ) const 
  { return Compare< MaximumFieldDimension >::not_equal( this->m_stride , rhs.m_stride ); }

  void print( std::ostream & os, const EntityRank & rank, const Part & part, FieldArrayRank field_rank ) const;

  private:
  EntityKey m_key ;
  size_type m_stride[ MaximumFieldDimension ];
};

typedef std::vector<FieldRestriction> FieldRestrictionVector;

std::string print_restriction( 
    const FieldRestriction & restr,
    const EntityRank & rank,
    const Part & part,
    FieldArrayRank field_rank
    );

} // namespace mesh
} // namespace stk

#endif // stk_mesh_baseImpl_FieldRestriction_hpp
