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

#ifndef phdmesh_Entity_hpp
#define phdmesh_Entity_hpp

//----------------------------------------------------------------------

#include <limits>
#include <iosfwd>

#include <util/Span.hpp>
#include <mesh/Kernel.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------
/** Connection between mesh entities with attributes.
 *
 *  Concept is a relation consisting of 
 *    { ( ( DomainEntity , RangeEntity ) , RelationIdentifier ) }
 *
 *  Design is for the DomainEntity to own a subset of the relation,
 *    DomainEntity -> { ( RangeEntity , RelationIdentifier ) }
 */

class Connect {
private:
  entity_key_type m_key ;
  Entity        * m_entity ;
public:

  ~Connect() {}

  Connect() : m_key(0), m_entity(NULL) {}

  Connect( const Connect & r ) : m_key(r.m_key), m_entity(r.m_entity) {}

  Connect & operator = ( const Connect & r )
    { m_key = r.m_key ; m_entity = r.m_entity ; return *this ; }

  /** Construct with range entity, identifier, and if anonymous */
  Connect( Entity & , unsigned );

  Entity * entity() const { return m_entity ; }

  unsigned       entity_type() const { return entity_rank( m_key ); }
  entity_id_type identifier()  const { return entity_id( m_key ); }

  entity_key_type key() const { return m_key ; }

  bool operator == ( const Connect & r ) const
    { return m_key == r.m_key && m_entity == r.m_entity ; }

  bool operator != ( const Connect & r ) const
    { return m_key != r.m_key || m_entity != r.m_entity ; }

  bool operator < ( const Connect & r ) const ;
};

typedef std::vector<Connect> ConnectSet ;

/** Span of a sorted connections for a given domain entity.
 *  Members are sorted by
 *  (1) range entity type,
 *  (2) connection type,
 *  (3) identifier, and
 *  (4) range entity identifier.
 */
typedef Span< ConnectSet::const_iterator > ConnectSpan ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** A mesh entity has an entity type, identifier, connections, and
 *  resides within a mesh kernel.  The mesh kernel holds its field data.
 */
class Entity : public SetvMember< entity_key_type > {
private:

  ConnectSet           m_connect ;    // Connections
  KernelSet::iterator  m_kernel ;     // Containing kernel
  unsigned             m_kernel_ord ; // Ordinal in the kernel
  unsigned             m_owner_rank ; // Parallel owner rank
  EntityProcSpan       m_sharing ;

public:

  /** Entity type */
  inline unsigned entity_type() const { return entity_rank( key() ); }

  /** Identifier */
  inline entity_id_type identifier() const { return entity_id( key() ); }

  /** Kernel in which this mesh entity resides */
  Kernel & kernel() const { return * m_kernel ; }

  /** Kernel in which this mesh entity resides */
  unsigned kernel_ordinal() const { return m_kernel_ord ; }

  /** Owning parallel rank */
  unsigned owner_rank() const { return m_owner_rank ; }

  //------------------------------------

  ConnectSpan connections() const { return ConnectSpan( m_connect ); }

  ConnectSpan connections( unsigned type ) const ;

  ConnectSpan connections( unsigned type , entity_id_type patch_id ) const ;

  //------------------------------------

  const EntityProcSpan & sharing() const { return m_sharing ; }

  //------------------------------------
  /** Pointer to field value for this mesh entity and field */
  template< class field_type >
  typename field_type::data_type * data( const field_type & f ) const
    { return m_kernel->data( f , m_kernel_ord ); }

  unsigned data_size( const FieldBase & f ) const
    { return m_kernel->data_size( f ); }

  template< class field_type >
  typename field_type::Dimension data_dim( const field_type & f ) const
    { return m_kernel->data_dim( f ); }

  //------------------------------------

  Entity();
  ~Entity();

private:

  Entity( const Entity & );
  Entity & operator = ( const Entity & );

  friend class Mesh ;
};

typedef Setv<Entity> EntitySet ;


//----------------------------------------------------------------------

/** Print identifier and connections */
std::ostream &
print_entity( std::ostream & , const std::string & lead , const Entity & );

std::ostream &
print_entity_key( std::ostream & , unsigned type , entity_id_type id );

std::ostream &
print_entity_key( std::ostream & , entity_key_type key );

std::ostream & print_connect( std::ostream & os ,
                              entity_key_type patch_key ,
                              entity_key_type entity_key );

std::ostream & operator << ( std::ostream & , const Connect & );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

