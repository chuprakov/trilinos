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

Mesh::Mesh( const Schema & schema ,
            ParallelMachine parallel ,
             unsigned kernel_capacity )
  : m_schema( schema ),
    m_parallel_machine( parallel ),
    m_parallel_size( parallel_machine_size( parallel ) ),
    m_parallel_rank( parallel_machine_rank( parallel ) ),
    m_kernel_capacity( kernel_capacity )
{
  static const char method[] = "phdmesh::Mesh::Mesh" ;

  m_schema.assert_committed( method );

  verify_parallel_consistency( schema , parallel );
}

Mesh::~Mesh()
{
  m_shares_all.clear();
  m_aura_domain.clear();
  m_aura_range.clear();

  // Remove entities from the kernels.
  // Destroy entities, which were allocated by the set itself.
  // Destroy kernels, which were *not* allocated by the set.

  for ( unsigned i = end_entity_rank ; 0 < i ; ) {
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
  for ( unsigned i = 0 ; i < end_entity_rank ; ++i ) {

    KernelSet & kset = m_kernels[ i ] ;

    for ( KernelSet::iterator ik = kset.begin() ; ik != kset.end() ; ++ik ) {
      ik->update_state();
    }
  }
}

//----------------------------------------------------------------------

const EntitySet & Mesh::entities( unsigned entity_type ) const
{
  const unsigned i = entity_type ;

  if ( end_entity_rank <= i ) {
    // Error
    std::string msg( "phdmesh::Mesh::entities FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_entities[ entity_type ];
}

const KernelSet & Mesh::kernels( unsigned entity_type ) const
{
  const unsigned i = entity_type ;

  if ( end_entity_rank <= i ) {
    // Error
    std::string msg( "phdmesh::Mesh::kernels FAILED with invalid type" );
    throw std::invalid_argument(msg);
  }

  return m_kernels[ entity_type ];
}

//----------------------------------------------------------------------

Entity * Mesh::get_entity( entity_key_type key ,
                           const char * required_by ) const
{
  const EntitySet & es = entities( entity_rank( key ) );
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

Entity & Mesh::declare_entity( entity_key_type key ,
                               const std::vector<Part*> & parts ,
                               int owner )
{
  const char method[] = "phdmesh::Mesh::declare_entity" ;

  if ( entity_id( key ) == 0 ) {
    std::string msg ;
    msg.append( method );
    msg.append( " FAILED : Not allowed to have identifier = zero" );
    throw std::logic_error( msg );
  }

  const unsigned entity_type = entity_rank( key );
  EntitySet   & eset = m_entities[ entity_type ] ;

  // Find or create the entity

  const std::pair< EntitySet::iterator , bool > result = eset.insert( key );

  const bool is_new = result.second ;

  if ( owner < 0 ) {
    owner = is_new ? parallel_rank() : result.first->owner_rank();
  }

  PartSet add( parts ) , remove ;

  if ( is_new && add.empty() ) {
    if ( owner == (int) parallel_rank() ) {
      Part * const owns_part = & m_schema.owns_part();
      add.push_back( owns_part );
    }
    else {
      Part * const univ_part = & m_schema.universal_part();
      add.push_back( univ_part );
    }
  }

  change_entity_parts( *result.first , add , remove );
  change_entity_owner( *result.first , (unsigned) owner );

  return * result.first ;
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

void Mesh::change_entity_identifier( Entity & e , entity_id_type id )
{
  static const char method[] = "phdmesh::Mesh::change_entity_identifier" ;

  const bool valid_id = id != 0 ;
  bool ok = valid_id ;

  if ( ok ) {
    const unsigned        type = e.entity_type();
    const entity_key_type key  = entity_key( type , id );
    Entity * const ptr = & e ;

    EntitySet & eset = m_entities[ type ] ;

    const std::pair< EntitySet::iterator , bool >
      result = eset.insert( key , ptr );

    ok = result.second || ( ptr == & * result.first );
  }

  if ( ! ok ) {
    std::ostringstream msg ;
    msg << "P" << parallel_size() ;
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

  if ( parallel_size() <= owner_rank ) {
    std::ostringstream msg ;
    msg << "P" << parallel_size() ;
    msg << ": " << method ;
    msg << "( " ;
    print_entity_key( msg , e.key() );
    msg << " , " ;
    msg << owner_rank ;
    msg << " ) FAILED due to invalid rank >= " ;
    msg << parallel_size() ;
    throw std::invalid_argument( msg.str() );
  }
  e.m_owner_rank = owner_rank ;
}

//----------------------------------------------------------------------

void Mesh::destroy_entity( Entity * e )
{
  while ( ! e->m_relation.empty() ) {
    destroy_relation( * e , * e->m_relation.back().entity() );
  }

  remove_entity( e->m_kernel , e->m_kernel_ord );

  e->m_kernel     = KernelSet::iterator();
  e->m_kernel_ord = 0 ;

  const unsigned  entity_type = e->entity_type();

  EntitySet & es = m_entities[ entity_type ];

  es.erase( e );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void verify_set_shares( const Mesh & M )
{
  static const char method[] = "phdmesh::Mesh::set_shares" ;

  std::string msg ;

  const EntityProcSet & shares = M.shares();

  // Parallel verification

  bool ok = comm_verify( M.parallel() , shares , msg );

  if ( ok ) {

    const unsigned p_rank = M.parallel_rank();
    Part & owns_part = M.schema().owns_part();
    Part & uses_part = M.schema().uses_part();

    std::ostringstream os ;

    os << "P" << p_rank << ": " << method << " FAILED " << std::endl ;

    const EntityProcSet::const_iterator es = shares.end();
          EntityProcSet::const_iterator is = shares.begin();

    for ( ; is != es ; ++is ) {

      // Verify not attempting to share with self

      if ( p_rank == is->second ) {
        print_entity_key( os , is->first->key() );
        os << " paired with P" << p_rank ;
        os << std::endl ;
        ok = false ;
      }

      // Verify each member has the uses part

      if ( ! is->first->kernel().has_superset( uses_part ) ) {
        print_entity_key( os , is->first->key() );
        os << " to be shared with P" ;
        os << is->second ;
        os << " does not have part " ;
        os << uses_part.name();
        os << std::endl ;
        ok = false ;
      }
    }

    // Verify all uses but not owned entities are shared
    // and their owner is one of the shared.

    for ( unsigned t = 0 ; t < end_entity_rank ; ++t ) {

      const KernelSet & kernels = M.kernels( t );

      const KernelSet::iterator ek = kernels.end();
            KernelSet::iterator ik = kernels.begin();

      for ( ; ek != ik ; ++ik ) {

        if ( ik->has_superset( uses_part ) ) {

          const Kernel::iterator e = ik->end();
                Kernel::iterator i = ik->begin();

          if ( ik->has_superset( owns_part ) ) {

            for ( ; e != i ; ++i ) {

              if ( (*i)->owner_rank() != p_rank ) {
                print_entity_key( os , (*i)->key() ); 
                os << " member of " ;
                os << owns_part.name();
                os << " but non-local owner = P" ;
                os << (*i)->owner_rank() ;
                os << std::endl ;
                ok = false ;
              }
            }
          }
          else {

            for ( ; e != i ; ++i ) {

              if ( (*i)->owner_rank() == p_rank ) {
                print_entity_key( os , (*i)->key() ); 
                os << " not member of " ;
                os << owns_part.name();
                os << " but local owner = P" ;
                os << p_rank ;
                os << std::endl ;
                ok = false ;
              }

              EntityProcSpan ss = (*i)->sharing();

              if ( ss.empty() ) {
                print_entity_key( os , (*i)->key() );
                os << " not member of " ;
                os << owns_part.name();
                os << ", is member of " ;
                os << uses_part.name();
                os << ", but is not shared." ;
                os << std::endl ;
                msg.append( os.str() );
                ok = false ;
              }

              for ( ; ss && (*i)->owner_rank() != ss->second ; ++ss );

              if ( ss.empty() ) {
                print_entity_key( os , (*i)->key() );
                os << " not member of " ;
                os << owns_part.name();
                os << ", is member of " ;
                os << uses_part.name();
                os << ", but does not share with owner P" ;
                os << (*i)->owner_rank();
                os << std::endl ;
                msg.append( os.str() );
                ok = false ;
              }
            }
          } // end not member of owns_spart
        }   // end not member of uses_part
      }     // end kernel loop
    }       // end entity type loop

    if ( ! ok ) { msg = os.str(); }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( M.parallel() , ReduceMin<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }
}

}

void Mesh::set_shares( const EntityProcSet & s )
{
  m_shares_all = s ;

  // Clear the entities' sharing in preparation for setting it.

  {
    const EntityProcSpan ss_empty ;

    for ( unsigned t = 0 ; t < end_entity_rank ; ++t ) {
      const EntitySet::iterator e = m_entities[t].end();
            EntitySet::iterator i = m_entities[t].begin();
      for ( ; e != i ; ++i ) { i->m_sharing = ss_empty ; }
    }
  }

  // Set the entities' sharing.

  const EntityProcSet::iterator es = m_shares_all.end();
        EntityProcSet::iterator is = m_shares_all.begin();

  while ( is != es ) {
    const EntityProcSet::iterator js = is ;
    for ( ; is != es && js->first == is->first ; ++is );
    js->first->m_sharing = EntityProcSpan( js , is );
  }

  // Verify 

  verify_set_shares( *this );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void verify_set_aura( const Mesh & M )
{
  static const char method[] = "phdmesh::Mesh::set_aura" ;

  std::string msg ;

  const EntityProcSet & domain = M.aura_domain();
  const EntityProcSet & range  = M.aura_range();

  // Parallel verification

  bool ok = comm_verify( M.parallel(), domain , range , msg );

  // Local verification of range ordering

  if ( ok ) { ok = verify( range , msg ); }

  if ( ok ) {

    Part & uses_part = M.schema().uses_part();
    Part & owns_part = M.schema().owns_part();

    const unsigned p_rank = M.parallel_rank();

    std::ostringstream os ;

    os << "P" << M.parallel_rank() << ": " << method << " FAILED " ;

    // Verify each entity in the domain has the owns part
    {
      Entity * e = NULL ;

      const EntityProcSet::const_iterator es = domain.end();
            EntityProcSet::const_iterator is = domain.begin();

      for ( ; is != es ; ++is ) {
        if ( e != is->first ) {
          e = is->first ;
        }

        if ( ! e->kernel().has_superset( owns_part ) ||
             p_rank != e->owner_rank() ) {
          print_entity_key( os , is->first->key() );
          os << " does not have part " ;
          os << owns_part.name();
          os << std::endl ;
          msg = os.str();
          ok = false ;
        }
      }
    }

    // Verify
    // (1) each entity in the range does not have the uses part
    // (2) an entity only appears in the range once

    {
      unsigned entity_count[ end_entity_rank ];
      unsigned kernel_count[ end_entity_rank ];

      for ( unsigned i = 0 ; i < end_entity_rank ; ++i ) {
        entity_count[i] = 0 ;
        kernel_count[i] = 0 ;
      }

      Entity * e = NULL ;

      const EntityProcSet::const_iterator es = range.end();
            EntityProcSet::const_iterator is = range.begin();

      for ( ; is != es ; ++is ) {
        const bool bad_domain = e == is->first ;
        e = is->first ;

        ++( entity_count[ e->entity_type() ] );

        const bool bad_part = e->kernel().has_superset( uses_part );
        const bool bad_rank = is->second != e->owner_rank();

        if ( bad_part || bad_rank || bad_domain ) {
          print_entity_key( os , e->key() );
          if ( bad_part ) {
            os << " member of " ;
            os << uses_part.name();
          }
          if ( bad_domain ) {
            os << " Received from multiple procesors " ;
          }
          if ( bad_rank ) {
            os << " Received from P" ;
            os << is->second ;
            os << " instead of owner P" ;
            os << e->owner_rank();
          }
          msg = os.str();
          ok = false ;
        }
      }

      // Verify aura range count equals not member of uses count.

      for ( unsigned i = 0 ; ok && i < end_entity_rank ; ++i ) {
        const KernelSet & kset = M.kernels( i );
        const KernelSet::const_iterator ek = kset.end();
              KernelSet::const_iterator ik = kset.begin();
        for ( ; ik != ek ; ++ik ) {
          if ( ! ik->has_superset( uses_part ) ) {
            kernel_count[i] += ik->size();
          }
        }

        if ( entity_count[i] != kernel_count[i] ) {
          os << " aura entity count = " << entity_count[i] ;
          os << " != " << kernel_count[i] << " = not 'uses' entity count" ;
          msg = os.str();
          ok = false ;
        }
      }
    }
  }

  {
    // Global reduce on the error flag, are any of the flag false ?
    unsigned flag = ok ;
    all_reduce( M.parallel() , ReduceMin<1>( & flag ) );
    ok = flag ;
  }

  if ( ! ok ) {
    throw std::runtime_error( msg );
  }

}

}

void Mesh::set_aura( const EntityProcSet & d , const EntityProcSet & r )
{
  m_aura_domain = d ;
  m_aura_range  = r ;

  verify_set_aura( *this );
}

//----------------------------------------------------------------------
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
  Mesh & mesh , Part & part , unsigned * const count )
{
  static const char method[] = "phdmesh::partset_entity_count" ;

  mesh.schema().assert_same_schema( method , part.schema() );

  for ( unsigned i = 0 ; i < end_entity_rank ; ++i ) {
    count[i] = 0 ;

    const KernelSet & ks = mesh.kernels( i );

    KernelSet::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( ik->has_superset( part ) ) {
        count[i] += ik->size();
      }
    }
  }
}

void partset_entity_count(
  Mesh & mesh , const PartSet & parts , unsigned * const count )
{
  static const char method[] = "phdmesh::partset_entity_count" ;

  for ( unsigned i = 0 ; i < end_entity_rank ; ++i ) {
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

    for ( unsigned i = 0 ; i < end_entity_rank ; ++i ) {
      const KernelSet & ks = mesh.kernels( i );

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

