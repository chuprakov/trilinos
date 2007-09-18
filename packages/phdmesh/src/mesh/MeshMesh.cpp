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

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Schema.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

Mesh::Mesh( const Schema & schema , const unsigned kernel_capacity[] )
  : m_schema( schema )
{
  static const char method[] = "phdmesh::Mesh::Mesh" ;

  m_schema.assert_committed( method );

  Copy<EntityTypeMaximum>( m_kernel_capacity , kernel_capacity );
}

Mesh::~Mesh()
{
  m_shares_all.clear();
  m_aura_domain.clear();
  m_aura_range.clear();

  // Remove entities from the kernels.
  // Destroy entities, which were allocated by the set itself.
  // Destroy kernels, which were *not* allocated by the set.

  for ( unsigned i = EntityTypeMaximum ; 0 < i ; ) {
    --i ;

    EntitySet & eset = m_entities[i] ;
    KernelSet & kset = m_kernels[i] ;

    for ( EntitySet::iterator ie = eset.begin() ; ie != eset.end() ; ++ie ) {
      ie->m_kernel = KernelSet::iterator();
    }
    eset.clear();

    while ( kset.size() ) { destroy_kernel( kset.begin() ); }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void Mesh::update_state()
{
  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {

    KernelSet & kset = m_kernels[ i ] ;

    for ( KernelSet::iterator ik = kset.begin() ; ik != kset.end() ; ++ik ) {
      ik->update_state();
    }
  }
}

//----------------------------------------------------------------------

const EntitySet & Mesh::entities( EntityType entity_type ) const
{
  const unsigned i = entity_type ;

  if ( EntityTypeMaximum <= i ) {
    // Error
    std::string msg( "phdmesh::Mesh::entities FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_entities[ entity_type ];
}

const KernelSet & Mesh::kernels( EntityType entity_type ) const
{
  const unsigned i = entity_type ;

  if ( EntityTypeMaximum <= i ) {
    // Error
    std::string msg( "phdmesh::Mesh::kernels FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_kernels[ entity_type ];
}

//----------------------------------------------------------------------

Entity * Mesh::get_entity( unsigned long key , const char * required_by ) const
{
  const EntitySet & es = entities( Entity::key_entity_type( key ) );
  EntitySet::iterator i = es.find( key );
  if ( required_by && i == es.end() ) {
    static const char method[] = "phdmesh::Mesh::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , key );
    msg << " , " << required_by << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
  return i != es.end() ? & * i : (Entity*) NULL ;
}

Entity * Mesh::get_entity( EntityType t , unsigned long id ,
                           const char * required_by ) const
{
  const EntitySet & es = entities( t );
  const unsigned long key = Entity::create_key( t , id );
  EntitySet::iterator i = es.find( key );
  if ( required_by && i == es.end() ) {
    static const char method[] = "phdmesh::Mesh::get_entity" ;
    std::ostringstream msg ;
    msg << method << "( " ;
    print_entity_key( msg , t , id );
    msg << " , " << required_by << " ) FAILED" ;
    throw std::runtime_error( msg.str() );
  }
  return i != es.end() ? & * i : (Entity*) NULL ;
}

Entity & Mesh::declare_entity( EntityType entity_type ,
                               unsigned long id ,
                               const std::vector<Part*> & parts ,
                               int owner )
{
  const char method[] = "phdmesh::Mesh::declare_entity" ;

  if ( id == 0 ) {
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED : Not allowed to have identifier = zero" );
    throw std::logic_error( msg );
  }

  unsigned long key  = Entity::create_key( entity_type , id );
  EntitySet   & eset = m_entities[ entity_type ] ;

  // Find or create the entity

  const std::pair< EntitySet::iterator , bool > result = eset.insert( key );

  {
    PartSet empty ;
    change_entity_parts( *result.first , parts , empty );
  }

  if ( result.second || 0 <= owner ) {
    if ( owner < 0 ) { owner = m_schema.parallel_rank(); }
    change_entity_owner( *result.first , (unsigned) owner );
  }

  return * result.first ;
}

Entity & Mesh::declare_entity( unsigned long key ,
                               const std::vector<Part*> & parts ,
                               int owner )
{
  const EntityType    t  = Entity::key_entity_type( key );
  const unsigned long id = Entity::key_identifier( key );
  return declare_entity( t , id , parts , owner );
}

//----------------------------------------------------------------------

void Mesh::change_entity_parts( Entity & e ,
                                const std::vector<Part*> & add_parts ,
                                const std::vector<Part*> & remove_parts )
{
  const char method[] = "phdmesh::Mesh::change_entity_parts" ;

  PartSet a_parts ;
  PartSet r_parts ;

  // If adding a part then add all supersets of the part

  for ( std::vector<Part*>::const_iterator
          ip =  add_parts.begin() ;
          ip != add_parts.end() ; ++ip ) {
    Part * const ap = *ip ;
    insert( a_parts , *ap );
 
    for ( std::vector<Part*>::const_iterator
            jp =  ap->supersets().begin() ;
            jp != ap->supersets().end() ; ++jp ) {
      insert( a_parts , **jp );
    }
  }

  // If removing a part then remove all subsets of that part.
  // This will include the 'intersection' parts.

  for ( std::vector<Part*>::const_iterator
          ip =  remove_parts.begin() ;
          ip != remove_parts.end() ; ++ip ) {
    Part * const rp = *ip ;
    insert( r_parts , *rp );
    for ( std::vector<Part*>::const_iterator
            jp =  rp->subsets().begin() ;
            jp != rp->subsets().end() ; ++jp ) {
      insert( r_parts , **jp );
    }
  }

  // Verify that r_parts and a_parts are disjoint and
  // that the remove parts does not contain the universal part.

  if ( intersect( a_parts , r_parts ) ||
       contain( r_parts , m_schema.universal_part() ) ) {
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED : Bad combination of add_parts and remove_parts" );
    msg.append( " add_parts = {" );
    for ( PartSet::iterator i = a_parts.begin() ; i != a_parts.end() ; ++i ) {
      msg.append( " " );
      msg.append( (*i)->name() );
    }
    msg.append( " }  remove_parts = {" );
    for ( PartSet::iterator i = r_parts.begin() ; i != r_parts.end() ; ++i ) {
      msg.append( " " );
      msg.append( (*i)->name() );
    }
    msg.append( " }" );
    throw std::logic_error( msg );
  }

  // Current kernel membership:

  bool changed = false ;

  std::vector<Part*> parts ;

  // Start with existing parts:

  const KernelSet::iterator ik_old = e.m_kernel ;

  if ( ik_old ) {

    ik_old->supersets( parts );

    // Remove parts:

    PartSet::iterator j = parts.begin();
    PartSet::iterator i = r_parts.begin() ;
    for ( ; i != r_parts.end() && j != parts.end() ; ++i ) {
      Part * const i_part = *i ;
      const unsigned i_ord = i_part->schema_ordinal();

      while ( j != parts.end() && (*j)->schema_ordinal() < i_ord ) { ++j ; }

      if ( j != parts.end() && *j == i_part ) {
        j = parts.erase( j );
        changed = true ;
      }
    }
  }

  // Add parts:
  for ( PartSet::iterator
          i = a_parts.begin() ;
          i != a_parts.end() ; ++i ) {
    Part & p = **i ;

    if ( ! contain( parts , p ) ) {
      // Add this part, check if any 'intersection' subsets
      // now fall within the collection of parts.

      insert( parts , p );

      for ( PartSet::const_iterator
              j = p.subsets().begin() ;
              j != p.subsets().end() ; ++j ) {
        Part & p_sub = **j ;
        if ( contain( parts , p_sub.intersection_of() ) ) {
          insert( parts , p_sub );
        }
      }

      changed = true ;
    }
  }

  if ( changed ) {
    insert_entity( e , parts );
  }
}

void Mesh::change_entity_identifier( Entity & e , unsigned long id )
{
  static const char method[] = "phdmesh::Mesh::change_entity_identifier" ;

  const bool valid_id = id != 0 ;
  bool ok = valid_id ;

  if ( ok ) {
    const EntityType    type = e.entity_type();
    const unsigned long key  = Entity::create_key( type , id );
    Entity * const ptr = & e ;

    EntitySet & eset = m_entities[ type ] ;

    const std::pair< EntitySet::iterator , bool >
      result = eset.insert( key , ptr );

    ok = result.second || ( ptr == & * result.first );
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "P" << m_schema.parallel_size() ;
    msg << ": " << method ;
    msg << "( " ;
    print_entity_key( msg , e.key() );
    msg << " , " ;
    msg << id ;
    msg << " ) FAILED" ;
    if ( valid_id ) { msg << " identifier is in use" ; }
    else            { msg << " identifier is invalid" ; }
    throw std::invalid_argument( msg.str() );
  }
}

void Mesh::change_entity_owner( Entity & e , unsigned owner_rank )
{
  static const char method[] = "phdmesh::Mesh::change_entity_owner" ;

  if ( m_schema.parallel_size() <= owner_rank ) {
    std::ostringstream msg ;
    msg << "P" << m_schema.parallel_size() ;
    msg << ": " << method ;
    msg << "( " ;
    print_entity_key( msg , e.key() );
    msg << " , " ;
    msg << owner_rank ;
    msg << " ) FAILED due to invalid rank >= " ;
    msg << m_schema.parallel_size() ;
    throw std::invalid_argument( msg.str() );
  }
  e.m_owner_rank = owner_rank ;
}

//----------------------------------------------------------------------

void Mesh::destroy_entity( Entity * e )
{
  const EntityType entity_type = e->entity_type();

  {
    const ConnectSet tmp( e->m_connect );
    const ConnectSet::const_iterator i_end = tmp.end();
          ConnectSet::const_iterator i     = tmp.begin();
    for ( ; i_end != i ; ++i ) {
      i->entity()->remove_connections( e );
    }
  }

  remove_entity( e->m_kernel , e->m_kernel_ord );

  e->m_kernel     = KernelSet::iterator();
  e->m_kernel_ord = 0 ;

  EntitySet & es = m_entities[ entity_type ];

  es.erase( e );
}

//----------------------------------------------------------------------

void Mesh::set_shares( const std::vector<EntityProc> & s )
{
  static const char method[] = "phdmesh::Mesh::set_shares" ;

  const unsigned p_rank = m_schema.parallel_rank();
  Part & shares_part = m_schema.shares_part();

  std::string msg ;

  // Parallel verification

  bool ok = comm_verify( m_schema.parallel() , s , msg );

  // Verify each member has the shares part

  unsigned entity_count[ EntityTypeMaximum ] ;
  unsigned kernel_count[ EntityTypeMaximum ] ;

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    entity_count[i] = 0 ;
    kernel_count[i] = 0 ;
  }
  Entity * e = NULL ;

  const std::vector<EntityProc>::const_iterator es = s.end();
        std::vector<EntityProc>::const_iterator is = s.begin();

  for ( ; ok && is != es ; ++is ) {
    if ( e != is->first ) {
      e = is->first ;
      ++( entity_count[ e->entity_type() ] );
    }

    if ( p_rank == is->second ) {
      std::ostringstream os ;
      os << "P" << p_rank << ": " << method << " FAILED " ;
      print_entity_key( os , is->first->key() );
      os << " paired with P" << p_rank ;
      msg = os.str();
      ok = false ;
    }

    if ( ! e->kernel().has_superset( shares_part ) ) {
      std::ostringstream os ;
      os << "P" << m_schema.parallel_rank() << ": " << method << " FAILED " ;
      print_entity_key( os , is->first->key() );
      os << " does not have part " ;
      os << shares_part.name();
      msg = os.str();
      ok = false ;
    }
  }

  if ( ok ) {
    // Verify that each entity with the shares part is in this array

    for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
      const KernelSet & kset = m_kernels[ EntityType(i) ];
      const KernelSet::const_iterator ek = kset.end();
            KernelSet::const_iterator ik = kset.begin();
      for ( ; ik != ek ; ++ik ) {
        if ( ik->has_superset( shares_part ) ) {
          kernel_count[i] += ik->size();
        }
      }
      if ( entity_count[i] != kernel_count[i] ) { ok = false ; }
    }

    if ( ! ok ) {
      std::ostringstream os ;
      os << "P" << m_schema.parallel_rank() << ": " << method ;
      os << " FAILED " << s.size();
      os << " entity counts {" ;
      for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
        os << " " << entity_count[i] ;
      }
      os << " } != kernel counts {" ;
      for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
        os << " " << kernel_count[i] ;
      }
      os << " }" ;
      msg = os.str();
    }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( m_schema.parallel() , ReduceMin<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }

  m_shares_all = s ;
}

//----------------------------------------------------------------------

void Mesh::set_aura( const std::vector<EntityProc> & d ,
                     const std::vector<EntityProc> & r )
{
  static const char method[] = "phdmesh::Mesh::set_aura" ;

  Part & owns_part = m_schema.owns_part();
  Part & aura_part = m_schema.aura_part();
  const unsigned p_rank = m_schema.parallel_rank();

  std::string msg ;

  // Parallel verification

  bool ok = comm_verify( m_schema.parallel() , d , r , msg );

  // Local verification of range ordering

  if ( ok ) { ok = verify( r , msg ); }

  // Verify each entity in the domain has the owns part
  {
    Entity * e = NULL ;

    const std::vector<EntityProc>::const_iterator es = d.end();
          std::vector<EntityProc>::const_iterator is = d.begin();

    for ( ; ok && is != es ; ++is ) {
      if ( e != is->first ) {
        e = is->first ;
      }

      if ( ! e->kernel().has_superset( owns_part ) ||
           p_rank != e->owner_rank() ) {
        std::ostringstream os ;
        os << "P" << m_schema.parallel_rank() << ": " << method << " FAILED " ;
        print_entity_key( os , is->first->key() );
        os << " does not have part " ;
        os << owns_part.name();
        msg = os.str();
        ok = false ;
      }
    }
  }

  // Verify each range member has the aura part
  {
    unsigned entity_count[ EntityTypeMaximum ];
    unsigned kernel_count[ EntityTypeMaximum ];

    for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
      entity_count[i] = 0 ;
      kernel_count[i] = 0 ;
    }

    Entity * e = NULL ;

    const std::vector<EntityProc>::const_iterator es = r.end();
          std::vector<EntityProc>::const_iterator is = r.begin();

    for ( ; ok && is != es ; ++is ) {
      const bool bad_domain = e == is->first ;
      e = is->first ;

      ++( entity_count[ e->entity_type() ] );

      const bool bad_part = ! e->kernel().has_superset( aura_part );
      const bool bad_rank = is->second != e->owner_rank();

      if ( bad_part || bad_rank || bad_domain ) {
        std::ostringstream os ;
        os << "P" << m_schema.parallel_rank() << ": " << method << " FAILED " ;
        print_entity_key( os , e->key() );
        if ( bad_domain ) {
          os << " Had multiple entries" ;
        }
        if ( bad_part ) {
          os << " Does not have part " ;
          os << aura_part.name();
        }
        if ( bad_rank ) {
          os << " Is received from P" ;
          os << is->second ;
          os << " instead of owner P" ;
          os << e->owner_rank();
        }
        msg = os.str();
        ok = false ;
      }
    }

    // Verify that each entity with the aura part is in this array

    for ( unsigned i = 0 ; ok && i < EntityTypeMaximum ; ++i ) {
      const KernelSet & kset = m_kernels[ EntityType(i) ];
      const KernelSet::const_iterator ek = kset.end();
            KernelSet::const_iterator ik = kset.begin();
      for ( ; ik != ek ; ++ik ) {
        if ( ik->has_superset( aura_part ) ) {
          kernel_count[i] += ik->size();
        }
      }

      if ( entity_count[i] != kernel_count[i] ) {
        std::ostringstream os ;
        os << "P" << m_schema.parallel_rank() << ": " << method << " FAILED " ;
        os << " aura entity count = " << entity_count[i] ;
        os << " != " << kernel_count[i] << " = aura kernel entity count" ;
        msg = os.str();
        ok = false ;
      }
    }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( m_schema.parallel() , ReduceMin<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }

  m_aura_domain = d ;
  m_aura_range = r ;
}

//----------------------------------------------------------------------

void get_kernels( const KernelSet & ks ,
                  Part & part , std::vector<Kernel*> & kernels )
{
  kernels.clear();

  const KernelSet::iterator ie = ks.end();
        KernelSet::iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    Kernel * const k = & * ik ;
    if ( k->has_superset( part ) ) {
      kernels.push_back( k );
    }
  }
}

void get_kernels_intersect(
  const KernelSet & ks ,
  const std::vector<Part*> & parts ,
  std::vector<Kernel*> & kernels )
{
  kernels.clear();

  const KernelSet::iterator ie = ks.end();
        KernelSet::iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    Kernel * const k = & * ik ;
    if ( k->has_superset( parts ) ) {
      kernels.push_back( k );
    }
  }
}

void get_kernels_union(
  const KernelSet & ks ,
  const std::vector<Part*> & parts ,
  std::vector<Kernel*> & kernels )
{
  kernels.clear();

  const KernelSet::iterator ie = ks.end();
        KernelSet::iterator ik = ks.begin();

  for ( ; ik != ie ; ++ik ) {
    Kernel * const k = & * ik ;
    bool result = false ;
    const std::vector<Part*>::const_iterator ep = parts.end();
          std::vector<Part*>::const_iterator ip = parts.begin();
    for ( ; ! result && ip != ep ; ++ip ) {
      result = k->has_superset( **ip );
    }
    if ( result ) {
      kernels.push_back( k );
    }
  }
}

//----------------------------------------------------------------------

void partset_entity_count(
  Mesh & mesh , Part & part , unsigned long * const count )
{
  static const char method[] = "phdmesh::partset_entity_count" ;

  mesh.schema().assert_same_schema( method , part.schema() );

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    count[i] = 0 ;

    const KernelSet & ks = mesh.kernels(  EntityType(i) );

    KernelSet::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( ik->has_superset( part ) ) {
        count[i] += ik->size();
      }
    }
  }
}

void partset_entity_count(
  Mesh & mesh , const PartSet & parts , unsigned long * const count )
{
  static const char method[] = "phdmesh::partset_entity_count" ;

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    count[i] = 0 ;
  }

  if ( ! parts.empty() ) {
    PartSet tmp( parts );

    order( tmp );

    {
      const PartSet::iterator j = tmp.end();
            PartSet::iterator i = tmp.begin();
      for ( ; i != j ; ++i ) {
        mesh.schema().assert_same_schema( method , (*i)->schema() );
      }
    }

    for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
      const KernelSet & ks = mesh.kernels(  EntityType(i) );

      KernelSet::const_iterator ik ;

      for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
        if ( ik->has_superset( tmp ) ) {
          count[i] += ik->size();
        }
      }
    }
  }
}

//----------------------------------------------------------------------

} // namespace phdmesh

