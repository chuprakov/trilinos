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

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>
#include <mesh/EntityManager.hpp>


namespace phdmesh {

//----------------------------------------------------------------------

void comm_mesh_remove_aura( Mesh & M )
{
  const Schema & S = M.schema();

  Part & aura_part = S.aura_part();

  unsigned k = EntityTypeMaximum ;

  while ( k ) {
    --k ;

    const KernelSet & ks = M.kernels( EntityType(k) );

    const KernelSet::const_iterator ek = ks.rend();
          KernelSet::const_iterator ik = ks.rbegin();

    while ( ik != ek ) {
      const Kernel & kernel = *ik ; --ik ;
      if ( kernel.has_superset( aura_part ) ) {
        unsigned ordinal ;
        do {
          ordinal = kernel.size() - 1 ;
          M.destroy_entity( kernel[ ordinal ] );
        } while ( ordinal );
      }
    }

    for ( ik = ks.rbegin() ; ik != ek ; ) {
      const Kernel & kernel = *ik ; --ik ;
      if ( kernel.has_superset( aura_part ) ) {
        std::string msg("phdmesh::comm_mesh_remove_aura impossible failure");
        throw std::runtime_error(msg);
      }
    }
  }

  const std::vector<EntityProc> tmp ;
  M.set_aura( tmp , tmp );
}

//----------------------------------------------------------------------
// Regenerate the parallel aura of mesh entities that are
// connected to the shared node mesh entities.

namespace {

bool not_member( const std::vector<EntityProc> & aura_domain ,
                 const EntityProc & ep )
{
  const std::vector<EntityProc>::const_iterator
    i = lower_bound( aura_domain , ep );

  return i == aura_domain.end() || ep != *i ;
}

class RegenAuraManager : public EntityManager {
private:
  RegenAuraManager( const RegenAuraManager & );
  RegenAuraManager & operator = ( const RegenAuraManager & );
public:
  RegenAuraManager() {}
  ~RegenAuraManager() {}

  const char * name() const ;

  void receive_entity(
    CommBuffer & buffer ,
    Mesh & receive_mesh ,
    const unsigned send_source ,
    std::vector<EntityProc> & receive_info ) const ;
};

const char * RegenAuraManager::name() const
{ static const char n[] = "phdmesh::RegenAuraManager" ; return n ; }

void RegenAuraManager::receive_entity(
  CommBuffer              & buffer ,
  Mesh                    & receive_mesh ,
  const unsigned            send_source ,
  std::vector<EntityProc> & receive_info ) const
{
  EntityType            entity_type ;
  unsigned long         entity_id ;
  unsigned              owner_rank ;
  std::vector<Part*>    parts ;
  std::vector<Connect>  connections ;
  std::vector<unsigned> send_destinations ;

  unpack_entity( buffer , receive_mesh ,
                 entity_type , entity_id , owner_rank ,
                 parts , connections , send_destinations );

  const Schema & S = receive_mesh.schema();
  Part * const owns_part   = & S.owns_part();
  Part * const shares_part = & S.shares_part();
  Part * const aura_part   = & S.aura_part();

  // Don't set the owns or shares part

  std::vector<Part*>::iterator ip = parts.end();
  while ( ip != parts.begin() ) {
    --ip ;
    Part * const p = *ip ;
    if ( p == owns_part || p == shares_part || p == aura_part ) {
      ip = parts.erase( ip );
    }
  }

  EntityProc ep ;

  ep.first = receive_mesh.get_entity( entity_type , entity_id );

  if ( NULL != ep.first ||
       send_source != owner_rank ||
       send_source == receive_mesh.parallel_rank() ) {
    std::string msg( "phdmesh::RegenAuraManager::receive_entity FAILED" );
    throw std::logic_error( msg );
  }

  parts.push_back( aura_part );

  ep.first = declare_entity( receive_mesh , entity_type , entity_id ,
                             owner_rank , parts , connections );
  ep.second = owner_rank ;

  receive_info.push_back( ep );
}

void scrub( Mesh & M , std::vector<EntityProc> & new_domain ,
                       std::vector<EntityProc> & new_range )
{
  static const char method[] = "phdmesh::comm_mesh_regenerate_aura::scrub" ;

  const Schema & S = M.schema();
  Part & shares_part = S.shares_part();
  Part & aura_part   = S.aura_part();
  const unsigned p_rank = M.parallel_rank();
  const unsigned p_size = M.parallel_size();

  std::vector<EntityProc>::iterator i ;

  for ( i = new_range.end() ; i != new_range.begin() ; ) {
    --i ;
    Entity & aura_entity = *i->first ;

    if ( ! aura_entity.kernel().has_superset( aura_part ) ) {
      i->first = NULL ;
    }
    else {

      bool destroy_it = true ;

      const ConnectSpan aura_con = aura_entity.connections();

      std::vector<Connect>::const_iterator ic = aura_con.first ;

      for ( ; destroy_it && ic != aura_con.second ; ++ic ) {
        Entity & e = * ic->entity();

        if ( ic->type() == Uses ) {
          if ( e.kernel().has_superset( shares_part ) ) {
            destroy_it = false ;
          }
        }
        else if ( ic->type() == UsedBy ) {
          if ( e.kernel().has_superset( aura_part ) ) {
            destroy_it = false ;
          }
          else {
            std::ostringstream msg ;
            msg << "P" << p_rank ;
            msg << " " << method ;
            msg << " FAILED Aura " ;
            print_entity_key( msg , aura_entity.key() );
            msg << " UsedBy non-aura " ;
            print_entity_key( msg , e.key() );
            throw std::runtime_error( msg.str() );
          }
        }
      }

      if ( destroy_it ) {
        M.destroy_entity( i->first );
        i->first = NULL ;
      }
    }
  }

  // Inform the owners of the destroyed aura entities

  CommAll all( M.parallel() );

  bool change = false ;

  for ( i = new_range.begin() ; i != new_range.end() ; ++i ) {
    const unsigned char flag = NULL == i->first ;
    if ( flag ) { change = true ; }
    all.send_buffer( i->second ).skip<unsigned char>(1);
  }

  change = all.allocate_buffers( p_size / 4 , false , change );

  if ( change ) {

    for ( i = new_range.begin() ; i != new_range.end() ; ++i ) {
      const unsigned char flag = NULL == i->first ;
      all.send_buffer( i->second ).pack<unsigned char>(flag);
    }

    all.communicate();

    for ( i = new_domain.begin() ; i != new_domain.end() ; ++i ) {
      unsigned char flag ;
      all.recv_buffer( i->second ).unpack<unsigned char>(flag);
      if ( flag ) { i->first = NULL ; }
    }

    for ( i = new_domain.end() ; i != new_domain.begin() ; ) {
      --i ;
      if ( NULL == i->first ) { i = new_domain.erase(i); }
    }

    for ( i = new_range.end() ; i != new_range.begin() ; ) {
      --i ;
      if ( NULL == i->first ) { i = new_range.erase(i); }
    }
  }
}

}


void comm_mesh_regenerate_aura( Mesh & M )
{
  static const char method[] = "phdmesh::comm_mesh_regenerate_aura" ;

  const RegenAuraManager mgr ;
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  const std::vector<EntityProc> & shares = M.shares();

  std::vector<EntityProc> old_aura_domain( M.aura_domain() );
  std::vector<EntityProc> old_aura_range(  M.aura_range() );

  scrub( M , old_aura_domain , old_aura_range );

  // Run the shared entities and push all 'UsedBy' entities
  // into the send vector.

  std::vector<EntityProc> new_aura_domain ;
  std::vector<EntityProc> new_aura_range ;

  for ( std::vector<EntityProc>::const_iterator
        i = shares.begin() ; i != shares.end() ; ) {

    const std::vector<EntityProc>::const_iterator ib = i ;
    for ( ; i != shares.end() && ib->first == i->first ; ++i );
    const std::vector<EntityProc>::const_iterator ie = i ;

    Entity * const entity = ib->first ;

    const ConnectSpan & entity_con = entity->connections();

    std::vector<Connect>::const_iterator ic ;

    for ( ic = entity_con.first ; ic != entity_con.second ; ++ic ) {

      Entity * const e = ic->entity();

      if ( ic->type() == UsedBy && p_rank == e->owner_rank() ) {

        const ConnectSpan e_con = e->connections();

        // Send to each shares processor

        std::vector<EntityProc>::const_iterator j ;

        for ( j = ib ; j != ie ; ++j ) {
          EntityProc entry ;

          entry.first = e ;
          entry.second = j->second ;

          // The entity->UsedBy->entity

          if ( not_member( shares , entry ) &&
               not_member( old_aura_domain , entry ) ) {
            new_aura_domain.push_back( entry );
          }

          // and all of the entity->UsedBy->entity->Uses->{entities}

          std::vector<Connect>::const_iterator jc ;

          for ( jc = e_con.first ; jc != e_con.second ; ++jc ) {
            if ( jc->type() == Uses && entity != jc->entity() ) {
              entry.first = jc->entity();
              if ( not_member( shares , entry ) ) {
                if ( p_rank == entry.first->owner_rank() ) {
                  if ( not_member( old_aura_domain , entry ) ) {
                    new_aura_domain.push_back( entry );
                  }
                }
                else if ( entry.second != entry.first->owner_rank() ) {
                  // Inform the owner of this need
                  new_aura_range.push_back( entry );
                }
              }
            }
          }
        }
      }
    }
  }

  {
    CommAll all( M.parallel() );

    std::vector<EntityProc>::iterator j ;

    for ( j = new_aura_range.begin() ; j != new_aura_range.end() ; ++j ) {
      all.send_buffer( j->first->owner_rank() ).skip<unsigned long>(2);
    }

    all.allocate_buffers( p_size / 4 , false /* not symmetric */ );

    for ( j = new_aura_range.begin() ; j != new_aura_range.end() ; ++j ) {
      unsigned long data[2] ;
      data[0] = j->first->key();
      data[1] = j->second ;
      all.send_buffer( j->first->owner_rank() ).pack<unsigned long>(data,2);
    }

    new_aura_range.clear();

    all.communicate();

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      CommBuffer & buf = all.recv_buffer(p);
      while ( buf.remaining() ) {
        unsigned long data[2] ;
        buf.unpack<unsigned long>( data , 2 );
        EntityProc tmp ;
        tmp.first  = M.get_entity( data[0] , method );
        tmp.second = data[1] ;
        if ( not_member( old_aura_domain , tmp ) ) {
          new_aura_domain.push_back( tmp );
        }
      }
    }
  }

  sort_unique( new_aura_domain );

  comm_mesh_entities( mgr , M , M , new_aura_domain , new_aura_range , false );

  new_aura_domain.insert( new_aura_domain.end() , old_aura_domain.begin() ,
                                                  old_aura_domain.end() );

  new_aura_range.insert( new_aura_range.end() , old_aura_range.begin() ,
                                                old_aura_range.end() );

  sort_unique( new_aura_domain );
  sort_unique( new_aura_range );

  M.set_aura( new_aura_domain , new_aura_range );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

