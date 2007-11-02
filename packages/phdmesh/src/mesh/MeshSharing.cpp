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
#include <util/ParallelIndex.hpp>
#include <util/ParallelReduce.hpp>

#include <mesh/Mesh.hpp>
#include <mesh/Schema.hpp>
#include <mesh/Comm.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_verify_shared_entity_values(
  const Mesh & M ,
  const Field<void,0> & f )
{
  const EntityType t = f.entity_type();
  const unsigned max_size   = f.max_size();
  const unsigned parallel_size = M.parallel_size();

  const std::pair< std::vector<EntityProc>::const_iterator ,
                   std::vector<EntityProc>::const_iterator >
    shares = span( M.shares() , t );

  std::vector<EntityProc>::const_iterator ic ;

  ParallelMachine comm = M.parallel();

  CommAll sparse ;

  {
    const unsigned zero = 0 ;
    std::vector<unsigned> comm_size( parallel_size , zero );

    for ( ic = shares.first ; ic != shares.second ; ++ic ) {
      comm_size[ ic->second ] += max_size ;
    }

    const unsigned * const p_size = & comm_size[0] ;

    sparse.allocate_buffers( comm, parallel_size / 4 , p_size, p_size );
  }

  std::vector<unsigned char> scratch( max_size );

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    Kernel & k = e.kernel();
    const FieldDimension & dim = dimension( f , k );
    const unsigned this_size = dim.size();
    CommBuffer & b = sparse.send_buffer( ic->second );

    unsigned char * data = (unsigned char *) e.data( f );

    b.pack<unsigned char>( data , this_size );
    b.skip<unsigned char>( max_size - this_size );
  }

  sparse.communicate();

  unsigned ok = 1 ;

  for ( ic = shares.first ; ic != shares.second ; ++ic ) {
    Entity & e = * ic->first ;
    Kernel & k = e.kernel();
    const FieldDimension & dim = dimension( f , k );
    const unsigned this_size = dim.size();
    CommBuffer & b = sparse.recv_buffer( ic->second );

    unsigned char * data = (unsigned char *) e.data( f );
    unsigned char * scr  = & scratch[0] ;

    b.unpack<unsigned char>( scr , this_size );
    b.skip<unsigned char>( max_size - this_size );

    // Compare data and scratch
    for ( unsigned j = 0 ; ok && j < this_size ; ++j ) {
      ok = data[j] == scr[j] ;
    }
  }

  all_reduce( comm , ReduceMin<1>( & ok ) );

  return ok ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void comm_mesh_discover_sharing( Mesh & M )
{
  static const char method[] = "phdmesh::comm_mesh_discover_sharing" ;

  const Schema & S = M.schema();
  const unsigned p_rank = M.parallel_rank();

  std::vector<unsigned long> local ;
  std::vector< ParallelIndex::KeyProc > global ;
  std::vector<EntityProc> share ;

  unsigned count = 0 ;
  for ( unsigned k = 0 ; k < EntityTypeMaximum ; ++k ) {
    const EntityType type = EntityType(k) ;
    const EntitySet & eset = M.entities( type );
    count += eset.size();
  }

  local.reserve( count );

  for ( unsigned k = 0 ; k < EntityTypeMaximum ; ++k ) {
    const EntityType type = EntityType(k) ;
    const EntitySet & eset = M.entities( type );

    const EntitySet::const_iterator e = eset.end();
          EntitySet::const_iterator i ;
    for ( i = eset.begin() ; i != e ; ++i ) { local.push_back( i->key() ); }
  }

  ParallelIndex par_index( M.parallel() , local );

  par_index.query( global );

  for ( std::vector< ParallelIndex::KeyProc >::iterator
        i = global.begin() ; i != global.end() ; ) {

    const unsigned long key = i->first ;

    EntityProc ep ;
    ep.first = M.get_entity( key , method );

    for ( ; i != global.end() && key == i->first ; ++i ) {
      ep.second = i->second ;
      share.push_back( ep );
    }
  }

  sort_unique( share );

  // Now revise the shared and ownership parts for the entities

  Part * const owns_part  = & S.owns_part();
  Part * const shares_part = & S.shares_part();

  const std::vector<EntityProc>::iterator ipe = share.end();
        std::vector<EntityProc>::iterator ip ;

  for ( ip = share.begin() ; ip != ipe ; ) {
    Entity * const entity = ip->first ;
    const unsigned p_send = ip->second ;

    for ( ; ip != ipe && ip->first == entity ; ++ip );

    std::vector<Part*> add_parts , remove_parts ;

    add_parts.push_back( shares_part );

    const unsigned p_owner = p_rank < p_send ? p_rank : p_send ;

    if ( p_owner != p_rank ) {
      remove_parts.push_back( owns_part );
    }
    else {
      add_parts.push_back( owns_part );
    }

    M.change_entity_parts( *entity , add_parts , remove_parts );
    M.change_entity_owner( *entity , p_owner );
  }

  // Set the shares connection
 
  M.set_shares( share );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// If an entity is not owned and not UsedBy an owned entity then remove it.

bool comm_mesh_scrub_sharing( Mesh & M )
{
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  std::vector<EntityProc> shares( M.shares() );

  bool changed = false ;

  std::vector<EntityProc>::iterator i ;

  for ( i = shares.end() ; i != shares.begin() ; ) {
    const std::vector<EntityProc>::iterator ie = i ;

    Entity * const entity = (--i)->first ;

    for ( ; i != shares.begin() && entity == i->first ; --i );
    if ( entity != i->first ) { ++i ; }

    bool destroy_it = p_rank != entity->owner_rank();

    const ConnectSpan con = entity->connections();
    std::vector<Connect>::const_iterator j = con.first ;
    for ( ; j != con.second && destroy_it ; ++j ) {
      destroy_it = ! ( UsedBy == j->type() &&
                       p_rank == j->entity()->owner_rank() );
    }

    if ( destroy_it ) {
      std::vector<EntityProc>::iterator ib = i ;
      for ( ; ib != ie ; ++ib ) { ib->first = NULL ; }
      M.destroy_entity( entity );
      changed = true ;
    }
  }

  // Inform sharing processors of the destruction

  CommAll all( M.parallel() );

  for ( i = shares.begin() ; i != shares.end() ; ++i ) {
    all.send_buffer( i->second ).skip<unsigned char>(1);
  }

  const bool symmetric = true ;

  changed = all.allocate_buffers( p_size / 4 , symmetric , changed );

  if ( changed ) {
    Part * const shares_part = & M.schema().shares_part();
    PartSet remove_parts ; remove_parts.push_back( shares_part );
    PartSet add_parts ;

    for ( i = shares.begin() ; i != shares.end() ; ++i ) {
      const unsigned char flag = NULL != i->first ;
      all.send_buffer( i->second ).pack<unsigned char>( flag );
    }

    all.communicate();

    for ( i = shares.begin() ; i != shares.end() ; ) {
      Entity * const e = i->first ;
      bool not_shared = true ;
      for ( ; i != shares.end() && e == i->first ; ++i ) {
        unsigned char flag ; 
        all.recv_buffer( i->second ).unpack<unsigned char>( flag );
        if ( flag ) { not_shared = false ; }
        else        { i->first = NULL ; }
      }
      if ( NULL != e && not_shared ) {
        // Remove the shared part
        M.change_entity_parts( *e , add_parts , remove_parts );
      }
    }
  
    for ( i = shares.end() ; i != shares.begin() ; ) {
      if ( NULL == (--i)->first ) { i = shares.erase( i ); }
    }

    M.set_shares( shares );

    comm_mesh_regenerate_aura( M );
  }

  return changed ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

void pack_info( CommAll & all , const std::vector<EntityProc> & shares )
{
  const std::vector<EntityProc>::const_iterator i_end = shares.end();
        std::vector<EntityProc>::const_iterator i     = shares.begin();

  while ( i != i_end ) {
    const EntityProc & ep = *i ; ++i ;

    Entity & entity = * ep.first ;
    const ConnectSpan connect = entity.connections();
    const unsigned long key = entity.key();
    const unsigned p_owner  = entity.owner_rank();
    const unsigned numconnect = std::distance( connect.first, connect.second );

    const unsigned p_send = ep.second ;

    std::vector<unsigned> part_ordinals ;
    entity.kernel().supersets( part_ordinals );
    const unsigned numparts = part_ordinals.size();

    CommBuffer & b = all.send_buffer( p_send );

    b.pack<unsigned long>( key );
    b.pack<unsigned>( p_owner );
    b.pack<unsigned>( numparts );
    b.pack<unsigned>( numconnect );

    b.pack<unsigned>( & part_ordinals[0] , numparts );

    for ( std::vector<Connect>::const_iterator k = connect.first ;
           k != connect.second ; ++k ) {
      unsigned long con[2] ;
      con[0] = k->attribute();
      con[1] = k->entity()->key();
      b.pack<unsigned long>( con , 2 );
    }
  }
}

bool unpack_info_verify(
  CommAll & all ,
  const std::vector<EntityProc> & shares ,
  std::string & error_msg )
{
  static const char method[] =
    "phdmesh::comm_mesh_verify_parallel_consistency" ;

  const unsigned u_zero = 0 ;
  const unsigned long ul_zero = 0 ;
  bool result = true ;

  const std::vector<EntityProc>::const_iterator i_end = shares.end();
        std::vector<EntityProc>::const_iterator i     = shares.begin();

  std::vector<unsigned> recv_ordinal ;
  std::vector<unsigned long> recv_connect ;

  while ( result && i != i_end ) {
    const EntityProc & ep = *i ; ++i ;

    Entity & entity = * ep.first ;
    Kernel & kernel = entity.kernel();
    Mesh   & mesh   = kernel.mesh();
    const Schema & schema = mesh.schema();
    const PartSet & mesh_parts = schema.get_parts();
    const ConnectSpan connect = entity.connections();
    Part * const owns_part  = & schema.owns_part();
    Part * const shares_part = & schema.shares_part();
    Part * const aura_part   = & schema.aura_part();

    const unsigned connect_size = std::distance(connect.first, connect.second);
    const unsigned long key = entity.key();
    const unsigned p_owner  = entity.owner_rank();
    const unsigned p_local  = mesh.parallel_rank();
    const unsigned p_recv   = ep.second ;

    std::vector<unsigned> part_ordinals ;
    kernel.supersets( part_ordinals );
    const unsigned parts_size = part_ordinals.size();

    ConnectSpan::first_type ic ;

    CommBuffer & b = all.recv_buffer( p_recv );

    unsigned long recv_key ;   b.unpack<unsigned long>( recv_key );
    unsigned recv_p_owner ;    b.unpack<unsigned>( recv_p_owner );
    unsigned recv_parts_size ; b.unpack<unsigned>( recv_parts_size );
    unsigned recv_connect_size ; b.unpack<unsigned>( recv_connect_size );

    recv_ordinal.assign( recv_parts_size , u_zero );
    {
      unsigned * const tmp = & recv_ordinal[0] ;
      b.unpack<unsigned>( tmp , recv_parts_size );
    }

    recv_connect.assign( recv_connect_size * 2 , ul_zero );
    {
      unsigned long * const tmp = & recv_connect[0] ;
      b.unpack<unsigned long>( tmp , recv_connect_size * 2 );
    }

    result = recv_key == key && recv_p_owner == p_owner ;

    if ( result ) {

      bool has_shared = true ;

      unsigned j_this = 0 ;
      unsigned j_recv = 0 ;

      while ( result &&
              j_this < parts_size &&
              j_recv < recv_parts_size ) {

        const int ord_recv = recv_ordinal[ j_recv ];
        if ( ord_recv < 0 || ((int) mesh_parts.size()) <= ord_recv ) {
          // Bad ordinal
          result = false ;
        }
        else {
          Part * const part_this = mesh_parts[ part_ordinals[ j_this ] ] ;
          Part * const part_recv = mesh_parts[ ord_recv ];

          if ( owns_part == part_this ) {
            // Local ownership by this processor
            result = p_local == p_owner &&
                     p_local == recv_p_owner ;
            ++j_this ;
          }
          else if ( owns_part == part_recv ) {
            // Remote ownership by sending processor
            result = p_recv == p_owner &&
                     p_recv == recv_p_owner ;
            ++j_recv ;
          }
          else if ( part_this != part_recv ) {
            result = false ;
          }
          else {
            ++j_this ;
            ++j_recv ;

            if ( shares_part == part_this ) { has_shared = true ; }
          }
        }
      }
      if ( ! has_shared ) { result = false ; }
    }

    if ( connect_size != recv_connect_size ) { result = false ; }

    ic = connect.first ;
    for ( unsigned k = 0 ; result && k < recv_connect_size ; ++k , ++ic ) {
      const Connect & con = *ic ;
      const unsigned k2 = k * 2 ;
      const unsigned      recv_con_attr = recv_connect[k2] ;
      const unsigned long recv_con_key  = recv_connect[k2+1] ;
      if ( recv_con_attr != con.attribute() ||
           recv_con_key  != con.entity()->key() ) { result = false ; }
    }

    if ( ! result ) {
      std::ostringstream os ;

      os << std::endl ;
      os << method << " ERROR :" << std::endl ;
      os << "  P" << p_local << " : " ;
      print_entity_key( os , key );
      os << ".{ Owner = P" << p_owner << " ; Parts =" ;
      for ( unsigned j = 0 ; j < part_ordinals.size() ; ++j ) {
        os << " " << mesh_parts[ part_ordinals[j] ]->name();
      }
      os << " ;" << std::endl << "  Connections =" ;
      for ( ic = connect.first ; ic != connect.second ; ++ic ) {
        Entity & con_e = * ic->entity();

        const bool has_owned  = con_e.kernel().has_superset( *owns_part );
        const bool has_shared = con_e.kernel().has_superset( *shares_part );
        const bool has_aura   = con_e.kernel().has_superset( *aura_part );

        os << std::endl << "    " ;
        os << *ic ;
        os << ".{ Owner = P" << con_e.owner_rank();
        if ( has_owned )  { os << " owned" ; }
        if ( has_shared ) { os << " shares" ; }
        if ( has_aura )   { os << " aura" ; }
        os << " }" ;
      }
      os << " } != " << std::endl ;
      os << "  P" << p_recv << " : " ;
      print_entity_key( os , recv_key );
      os << ".{ Owner = P" << recv_p_owner << " ; Parts =" ;
      for ( unsigned j = 0 ; j < recv_parts_size ; ++j ) {
        const unsigned ord_recv = recv_ordinal[j] ;
        if ( mesh_parts.size() <= ord_recv ) {
          os << " unsigned[" << ord_recv << "]" ;
        }
        else {
          os << " " << mesh_parts[ord_recv]->name();
        }
      }
      os << " ;" << std::endl << "  Connections =" ;
      for ( unsigned j = 0 ; j < recv_connect_size * 2 ; ) {
        const unsigned      attr    = recv_connect[j] ; ++j ;
        const unsigned long con_key = recv_connect[j] ; ++j ;
        const Connect con( attr );
        os << std::endl << "    " ;
        os << connect_type_name( con.type() );
        os << "[" << con.identifier() << "]->" ;
        print_entity_key( os , con_key );
      }
      os << " }" << std::endl ;

      error_msg = os.str();
    }
  }

  return result ;
}

bool verify_parallel_attributes(
  Mesh & M , EntityType type , std::string & msg )
{
  bool result = true ;

  // Verify all entities with the shares_part are in the sharing list

  const Schema & S = M.schema();
  const unsigned p_rank = M.parallel_rank();
  Part & owns_part  = S.owns_part();
  Part & shares_part = S.shares_part();
  Part & aura_part   = S.aura_part();

  const EntitySet & es = M.entities( type );

  std::pair< std::vector<EntityProc>::const_iterator ,
             std::vector<EntityProc>::const_iterator >
    share_span = span( M.shares() , type ) ,
    aura_span = span( M.aura_range() , type );

  // Iterate everything in identifier ordering.

  const EntitySet::iterator i_end = es.end();
        EntitySet::iterator i     = es.begin();

  while ( result && i != i_end ) {
    Entity & entity = *i ; ++i ;
    Kernel & kernel = entity.kernel();

    const bool owned  = kernel.has_superset( owns_part );
    const bool shared = kernel.has_superset( shares_part );
    const bool aura   = kernel.has_superset( aura_part );
    const unsigned p_owner = entity.owner_rank();

    if ( owned ) {
      // Owner is local, cannot be aura
      if ( p_owner != p_rank ) { result = false ; }
      if ( aura ) { result = false ; }
    }
    else {
      // Owner is remote, must be shared or aura
      if ( p_owner == p_rank ) { result = false ; }
      if ( ! ( shared || aura ) ) { result = false ; }
    }

    // Cannot be shared and aura

    if ( aura && shared ) { result = false ; }

    bool in_shared = false ;
    bool in_aura = false ;
    unsigned aura_rank = p_owner ;

    // Could appear multiple times in the shared parallel connection:

    while ( share_span.first != share_span.second &&
            share_span.first->first == & entity ) {
      in_shared = true ;
      ++share_span.first ;
    }

    // Can appear at most once in the aura parallel connection:

    if ( aura_span.first != aura_span.second &&
         aura_span.first->first == & entity ) {

      in_aura = true ;

      aura_rank = aura_span.first->second ;

      if ( p_owner != aura_rank ) { result = false ; }

      ++aura_span.first ;
    }

    if ( shared != in_shared ) { result = false ; }
    if ( aura   != in_aura )   { result = false ; }

    if ( ! result ) {
      std::ostringstream os ;
      os << "P" << p_rank << " : Error with parallel attributes "
         << std::endl ;
      os << "  " ;
      print_entity_key( os , entity.key() );
      os << ".owner(P" << p_owner << ")" ;
      if ( owned )  { os << " , owns" ; }
      if ( shared ) { os << " , shares" ; }
      if ( aura )   { os << " , aura" ; }
      if ( in_shared ) { os << " , in_shared" ; }
      else             { os << " , not_in_shared" ; }
      if ( in_aura ) { os << " , in_aura(P" << aura_rank << ")" ; }
      else           { os << " , not_in_aura" ; }
      os << std::endl ;
      msg = os.str();
    }
  }
  return result ;
}

}

bool comm_mesh_verify_parallel_consistency( Mesh & M )
{
  // Exchange the shared entities' identifiers, part mesh ordinals, and
  // owner processor.  Should be fully consistent modulo the owner part.

  const unsigned p_size = M.parallel_size();

  std::string msg ;
  int result = 1 ;

  if ( result ) { result = verify_parallel_attributes( M , Node , msg ); }
  if ( result ) { result = verify_parallel_attributes( M , Edge , msg ); }
  if ( result ) { result = verify_parallel_attributes( M , Face , msg ); }
  if ( result ) { result = verify_parallel_attributes( M , Element , msg ); }
  if ( result ) { result = verify_parallel_attributes( M , Other , msg ); }

  {
    // Verify all shared entities

    const std::vector<EntityProc> & shares = M.shares();

    CommAll all_info( M.parallel() );

    pack_info( all_info , shares );

    all_info.allocate_buffers( p_size / 4 , false );

    pack_info( all_info , shares );

    all_info.communicate();

    if ( ! unpack_info_verify( all_info , shares , msg ) ) { result = false ; }
  }

  // Verify consistency of aura

  if ( ! comm_verify( M.parallel(), M.aura_domain(), M.aura_range(), msg) ) {
    result = false ;
  }

  // Global reduction of result flag:

  all_reduce( M.parallel() , ReduceMin<1>( & result ) );

  //--------------------------------------------------------------------
  // If an error occured collect the messages on
  // processor zero and output to standard error stream.

  if ( ! result ) { all_write( M.parallel() , std::cerr , msg ); }

  return result ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

class SharingManager : public EntityManager {
private:
  SharingManager( const SharingManager & );
  SharingManager & operator = ( const SharingManager & );
public:
  SharingManager() {}
  ~SharingManager() {}

  const char * name() const ;

  void receive_entity(
    CommBuffer & buffer ,
    Mesh & receive_mesh ,
    const unsigned send_source ,
    std::vector<EntityProc> & receive_info ) const ;
};

const char * SharingManager::name() const
{ static const char n[] = "phdmesh::SharingManager" ; return n ; }

void SharingManager::receive_entity(
  CommBuffer & buffer ,
  Mesh & receive_mesh ,
  const unsigned send_source ,
  std::vector<EntityProc> & receive_info ) const
{
  // receive_info is the new_aura_range

  const Schema & S = receive_mesh.schema();
  Part * const owns_part  = & S.owns_part();
  Part * const shares_part = & S.shares_part();
  Part * const aura_part   = & S.aura_part();

  EntityType            entity_type ;
  unsigned long         entity_id ;
  unsigned              owner_rank ;
  std::vector<Part*>    parts ;
  std::vector<Connect>  connections ;
  std::vector<unsigned> send_dest ;

  unpack_entity( buffer , receive_mesh ,
                 entity_type , entity_id , owner_rank ,
                 parts , connections , send_dest );

  if ( send_source != owner_rank ) {
    std::string msg( "phdmesh::SharingManager::receive_entity FAILED owner" );
    throw std::logic_error( msg );
  }

  { // Remove owned part, add shared part.
    std::vector<Part*>::iterator ip = parts.begin() ;
    for ( ; ip != parts.end() && owns_part != *ip ; ++ip );
    if ( ip != parts.end() ) { parts.erase( ip ); }
    parts.push_back( shares_part );
  }

  EntityProc ep ;

  ep.first  = receive_mesh.get_entity( entity_type , entity_id );
  ep.second = send_source ;

  // If entity exists it must be an aura

  if ( NULL != ep.first ) {
    std::vector<EntityProc>::iterator
      k = lower_bound( receive_info , *ep.first );

    if ( ! ep.first->kernel().has_superset( *aura_part ) ||
         k == receive_info.end() || k->first != ep.first ) {
      std::string msg( "phdmesh::SharingManager::receive_entity FAILED aura" );
      throw std::logic_error( msg );
    }

    receive_info.erase( k );

    std::vector<Part*> p_remove ; p_remove.push_back( aura_part );
    receive_mesh.change_entity_parts( *ep.first , parts , p_remove );
  }

  ep.first = declare_entity( receive_mesh , entity_type , entity_id ,
                             owner_rank , parts , connections );
}

bool not_member( const std::vector<EntityProc> & s , const EntityProc & m )
{
  const std::vector<EntityProc>::const_iterator i = lower_bound( s , m );
  return i == s.end() || m != *i ;
}

// Owner packs sharing

void pack_sharing( CommAll & all , const std::vector<EntityProc> & sharing )
{
  const unsigned p_rank = all.parallel_rank();

  std::vector<EntityProc>::const_iterator i , j ;

  for ( i = sharing.begin() ; i != sharing.end() ; ) {
    const std::vector<EntityProc>::const_iterator ib = i ;
    for ( ; i != sharing.end() && i->first == ib->first ; ++i );
    const std::vector<EntityProc>::const_iterator ie = i ;

    if ( p_rank == ib->first->owner_rank() ) {
      const unsigned n = std::distance( ib , ie ) - 1 ;
      const unsigned long key = ib->first->key();
      for ( i = ib ; i != ie ; ++i ) {
        CommBuffer & buf = all.send_buffer( i->second );
        buf.pack<unsigned long>( key );
        buf.pack<unsigned>( n );
        if ( n ) {
          for ( j = ib ; j != ie ; ++j ) {
            if ( i != j ) { buf.pack<unsigned>( j->second ); }
          }
        }
      }
    }
  }
}

void unpack_sharing( CommAll & all , Mesh & M ,
                     std::vector<EntityProc> & sharing ,
                     const char * method )
{
  const unsigned p_size = all.parallel_size();

  for ( unsigned p = 0 ; p < p_size ; ++p ) {
    CommBuffer & buf = all.recv_buffer( p );
    while ( buf.remaining() ) {
      unsigned long key ; buf.unpack<unsigned long>( key );
      unsigned n ;        buf.unpack<unsigned>( n );
      EntityProc tmp ;
      tmp.first = M.get_entity( key , method );
      tmp.second = p ; // The owner sent it
      sharing.push_back( tmp );
      for ( unsigned i = 0 ; i < n ; ++i ) {
        buf.unpack<unsigned>( tmp.second );
        sharing.push_back( tmp );
      }
    }
  }
}


}

void comm_mesh_add_sharing( Mesh & M , const std::vector<EntityProc> & add )
{
  static const char method[] = "phdmesh::comm_mesh_add_sharing" ;

  const SharingManager mgr ;
  const Schema & S = M.schema();
  const unsigned p_rank = M.parallel_rank();
  const unsigned p_size = M.parallel_size();
  Part * const shares_part = & S.shares_part();

  const std::vector<EntityProc> & old_shares = M.shares();

  std::vector<EntityProc> new_shares ;

  //--------------------------------------------------------------------
  // Have the owners send the to-be-shared entities to the
  // specified processors.  If an aura relationship exists
  // then remove it and replace it with a shared relationship.
  {
    std::vector<EntityProc> new_aura_domain( M.aura_domain() );
    std::vector<EntityProc> new_aura_range(  M.aura_range() );

    {
      CommAll all( M.parallel() );

      for ( std::vector<EntityProc>::const_iterator
            j = add.begin() ; j != add.end() ; ++j ) {

        const unsigned p_owner = j->first->owner_rank();

        if ( p_owner != j->second ) {
          if ( p_owner != p_rank ) {
            all.send_buffer( p_owner ).skip<unsigned long>(2);
          }
          else if ( not_member( old_shares , *j ) ) { 
            new_shares.push_back( *j );
          }
        }
      }

      all.allocate_buffers( p_size / 4 , false /* not symmetric */ );

      for ( std::vector<EntityProc>::const_iterator
            j = add.begin() ; j != add.end() ; ++j ) {

        const unsigned p_owner = j->first->owner_rank();

        if ( p_owner != j->second && p_owner != p_rank ) {
          unsigned long data[2] ;
          data[0] = j->first->key();
          data[1] = j->second ;
          all.send_buffer( p_owner ).pack<unsigned long>(data,2);
        }
      }

      all.communicate();

      for ( unsigned j = 0 ; j < p_size ; ++j ) {
        CommBuffer & buf = all.recv_buffer(j);
        while ( buf.remaining() ) {
          unsigned long data[2] ;
          buf.unpack<unsigned long>( data , 2 );
          EntityProc tmp ;
          tmp.first = M.get_entity( data[0] , method );
          tmp.second = data[1] ;
          if ( not_member( old_shares , tmp ) ) {
            new_shares.push_back( tmp );
          }
        }
      }
    }

    sort_unique( new_shares );

    // Change the newly shared mesh entities on the owner processor
    // Remove from the aura_domain.
    {
      Entity * e = NULL ;
      PartSet remove_parts ;
      PartSet add_parts ; add_parts.push_back( shares_part );
      std::vector<EntityProc>::iterator i ;
      for ( i = new_shares.begin() ; i != new_shares.end() ; ++i ) {
        if ( e != i->first ) {
          e = i->first ;
          M.change_entity_parts( *e , add_parts , remove_parts );
        }

        const std::vector<EntityProc>::iterator
          k = lower_bound( new_aura_domain , *i );
        if ( k != new_aura_domain.end() && *k == *i ) {
          new_aura_domain.erase( k );
        }
      }
    }

    // Communicate new-shared mesh entities from owner to sharing processors.
    // Remove from the aura_range, if present.

    comm_mesh_entities( mgr, M, M, new_shares, new_aura_range, false );

    // Some aura entities have become shared entities.

    M.set_aura( new_aura_domain , new_aura_range );
  }
  // Shared mesh entities are now in place.
  //--------------------------------------------------------------------
  // Owners inform the sharing processors of the extent of the sharing.

  new_shares.insert( new_shares.end() , old_shares.begin() ,
                                        old_shares.end() );

  sort_unique( new_shares );

  // Now owners send to sharing the extent of sharing
  {
    CommAll all( M.parallel() );

    pack_sharing( all , new_shares );

    all.allocate_buffers( p_size / 4 , false );

    pack_sharing( all , new_shares );

    all.communicate();

    unpack_sharing( all , M , new_shares , method );
  }

  sort_unique( new_shares );

  M.set_shares( new_shares );

  // Regenerate aura

  comm_mesh_regenerate_aura( M );
}


}

