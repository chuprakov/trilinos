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

#include <mesh/Kernel.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------
/** Intrinsic type of connections between mesh entities. */

enum ConnectType { Uses = 1 , UsedBy = 2 , Anonymous = 3 };

/** Query text name for entity connections type */

const char * connect_type_name( ConnectType );

//----------------------------------------------------------------------
/** Connection between mesh entities with attributes.
 *
 *  Concept is a relation consisting of 
 *    { ( ( DomainEntity , RangeEntity ) , Attribute ) }
 *
 *  Design is for the DomainEntity to own a subset of the relation,
 *    DomainEntity -> { ( RangeEntity , Attribute ) }
 *
 *  The Connect class consists of ( RangeEntity , Attribute ) where
 *  Attribute is the aggregation of ( EntityType , ConnectType , identifier ).
 */

class Connect {
private:
  unsigned m_attr ;
  Entity * m_entity ;

  enum { u_digits = std::numeric_limits<unsigned>::digits };
  enum { e_digits = EntityTypeDigits };
  enum { t_digits = 2 };
  enum { e_shift = u_digits - e_digits };
  enum { t_shift = e_shift  - t_digits };
  enum { e_mask = EntityTypeMask };
  enum { t_mask = 0x03 };
  enum { i_mask = ( ~((unsigned)0) ) >> ( e_digits + t_digits ) };
public:

  ~Connect() {}

  Connect() : m_attr(0), m_entity(NULL) {}

  Connect( const Connect & r ) : m_attr(r.m_attr), m_entity(r.m_entity) {}

  Connect & operator = ( const Connect & r )
    { m_attr = r.m_attr ; m_entity = r.m_entity ; return *this ; }

  /** Construct with range entity, connection type, and identifier */
  Connect( Entity & , ConnectType , unsigned );

  Entity * entity() const { return m_entity ; }

  ConnectType type() const
    { return ConnectType( ( m_attr >> t_shift ) & t_mask ); }

  unsigned identifier() const { return m_attr & i_mask ; }

  EntityType entity_type() const
    { return EntityType( ( m_attr >> e_shift ) & e_mask ); }

  bool operator == ( const Connect & r ) const
    { return m_attr == r.m_attr && m_entity == r.m_entity ; }

  bool operator != ( const Connect & r ) const
    { return m_attr != r.m_attr || m_entity != r.m_entity ; }

  bool operator < ( const Connect & r ) const ;

  //----------------------------------------
  // Methods for manipulation of the connection attribute.
  // Intended for internal use only.

  unsigned attribute() const { return m_attr ; }

  Connect( Entity & arg_entity , unsigned arg_attr )
    : m_attr( arg_attr ) , m_entity( & arg_entity ) {}

  explicit Connect( unsigned arg_attr )
    : m_attr( arg_attr ) , m_entity( NULL ) {}

  static unsigned attribute( unsigned entity_type ,
                             unsigned connect_type ,
                             unsigned identifier )
    {
      return ( ( entity_type  & e_mask ) << e_shift ) |
             ( ( connect_type & t_mask ) << t_shift ) |
               ( identifier   & i_mask ) ;
    }
};

typedef std::vector<Connect> ConnectSet ;

/** Span of a sorted connections for a given domain entity.
 *  Members are sorted by
 *  (1) range entity type,
 *  (2) connection type,
 *  (3) identifier, and
 *  (4) range entity identifier.
 */
typedef std::pair< ConnectSet::const_iterator ,
                   ConnectSet::const_iterator > ConnectSpan ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** A mesh entity has an entity type, identifier, connections, and
 *  resides within a mesh kernel.  The mesh kernel holds its field data.
 */
class Entity : public SetvMember<unsigned long> {
private:

  ConnectSet           m_connect ;    // Connections
  KernelSet::iterator  m_kernel ;     // Containing kernel
  unsigned             m_kernel_ord ; // Ordinal in the kernel
  unsigned             m_owner_rank ; // Parallel owner rank

  enum { KeyDigits = std::numeric_limits<unsigned long>::digits };
  enum { KeyIdentifierDigits = KeyDigits - EntityTypeDigits };
  enum { KeyIdentifierMask = ( ~((unsigned long)0) ) >> EntityTypeDigits };

public:

  inline static EntityType key_entity_type( unsigned long key )
    { return EntityType( key >> KeyIdentifierDigits ); }

  inline static unsigned long key_identifier( unsigned long key )
    { return key & KeyIdentifierMask ; }

  /** Entity type */
  inline EntityType entity_type() const { return key_entity_type(key()); }

  /** Identifier */
  inline unsigned long identifier() const { return key_identifier(key()); }

  /** Kernel in which this mesh entity resides */
  Kernel & kernel() const { return * m_kernel ; }

  /** Kernel in which this mesh entity resides */
  unsigned kernel_ordinal() const { return m_kernel_ord ; }

  /** Owning parallel rank */
  unsigned owner_rank() const { return m_owner_rank ; }

  //------------------------------------

  ConnectSpan connections() const
    { return ConnectSpan( m_connect.begin() , m_connect.end() ); }

  ConnectSpan connections( EntityType ) const ;

  ConnectSpan connections( EntityType , ConnectType ) const ;

  //------------------------------------
  /** Pointer to field value for this mesh entity and field */
  template<typename T,unsigned NDim> T * data( const Field<T,NDim> & f ) const
    { return m_kernel->data( f , m_kernel_ord ); }

  unsigned data_size( const Field<void,0> & f ) const
    { return m_kernel->data_size( f ); }

  //------------------------------------

  Entity();
  ~Entity();

  static unsigned long create_key( EntityType type , unsigned long id );

private:

  Entity( const Entity & );
  Entity & operator = ( const Entity & );

  friend class Mesh ;
};

typedef Setv<Entity> EntitySet ;

inline
const FieldDimension & dimension( const Field<void,0> & f , const Entity & e )
{ return dimension( f , e.kernel() ); }

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & , const Connect & );

/** Print identifier and connections */
std::ostream &
print_entity( std::ostream & , const std::string & lead , const Entity & );

std::ostream &
print_entity_key( std::ostream & , EntityType type , unsigned long id );

std::ostream &
print_entity_key( std::ostream & , unsigned long key );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

