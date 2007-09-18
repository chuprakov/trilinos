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
#include <sstream>
#include <algorithm>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>
#include <mesh/EntityManager.hpp>


namespace phdmesh {

//----------------------------------------------------------------------

namespace {

void print_entry( std::ostream & msg , const EntityProc & e )
{
  msg << "( " ;
  if ( e.first == NULL ) {
     msg << "NULL" ;
  }
  else {
    print_entity_key( msg , e.first->key() );
  }
  msg << " , " << e.second << " )" ;
}

}

//----------------------------------------------------------------------

struct LessEntityProc {
  LessEntityProc() {}

  bool operator()( const EntityProc & lhs , const EntityProc & rhs ) const
    {
      const unsigned long lhs_key = lhs.first->key();
      const unsigned long rhs_key = rhs.first->key();
      return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
    }

  bool operator()( const EntityProc & lhs , const Entity & rhs ) const
    {
      const unsigned long lhs_key = lhs.first->key();
      const unsigned long rhs_key = rhs.key();
      return lhs_key < rhs_key ;
    }

  bool operator()( const EntityProc & lhs , const unsigned rhs ) const
    {
      const unsigned lhs_type = lhs.first->entity_type();
      return lhs_type < rhs ;
    }
};

struct EqualEntityProc {
  bool operator()( const EntityProc & lhs , const EntityProc & rhs ) const
    { return lhs.first == rhs.first && lhs.second == rhs.second ; }
};

//----------------------------------------------------------------------

void sort_unique( std::vector<EntityProc> & v )
{
  std::vector<EntityProc>::iterator i = v.begin();
  std::vector<EntityProc>::iterator e = v.end();

  std::sort( i , e , LessEntityProc() );
  i = std::unique( i , e , EqualEntityProc() );
  v.erase( i , e );
}

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & v , Entity & m )
{
  const std::vector<EntityProc>::const_iterator i = v.begin();
  const std::vector<EntityProc>::const_iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::const_iterator
lower_bound( const std::vector<EntityProc> & v , const EntityProc & m )
{
  const std::vector<EntityProc>::const_iterator i = v.begin();
  const std::vector<EntityProc>::const_iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & v , Entity & m )
{
  const std::vector<EntityProc>::iterator i = v.begin();
  const std::vector<EntityProc>::iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::vector<EntityProc>::iterator
lower_bound( std::vector<EntityProc> & v , const EntityProc & m )
{
  const std::vector<EntityProc>::iterator i = v.begin();
  const std::vector<EntityProc>::iterator e = v.end();
  return lower_bound( i , e , m , LessEntityProc() );
}

std::pair< std::vector<EntityProc>::const_iterator ,
           std::vector<EntityProc>::const_iterator >
span( const std::vector<EntityProc> & v , EntityType entity_type )
{
  const unsigned t1 = entity_type ;
  const unsigned t2 = t1 + 1 ;

  std::pair< std::vector<EntityProc>::const_iterator ,
             std::vector<EntityProc>::const_iterator > result ;

  result.first  = v.begin();
  result.second = v.end();
  result.first  = lower_bound(result.first,result.second,t1,LessEntityProc());
  result.second = lower_bound(result.first,result.second,t2,LessEntityProc());
  return result ;
}

std::pair< std::vector<EntityProc>::iterator ,
           std::vector<EntityProc>::iterator >
span( std::vector<EntityProc> & v , EntityType entity_type )
{
  const unsigned t1 = entity_type ;
  const unsigned t2 = t1 + 1 ;

  std::pair< std::vector<EntityProc>::iterator ,
             std::vector<EntityProc>::iterator > result ;

  result.first  = v.begin();
  result.second = v.end();
  result.first  = lower_bound(result.first,result.second,t1,LessEntityProc());
  result.second = lower_bound(result.first,result.second,t2,LessEntityProc());
  return result ;
}

//----------------------------------------------------------------------

void procs( const std::vector<EntityProc> & v , Entity & m ,
            std::vector<unsigned> & p )
{
  p.clear();
  std::vector<EntityProc>::const_iterator i = lower_bound( v , m );
  for ( ; i != v.end() && i->first == & m ; ++i ) {
    p.push_back( i->second );
  }
}

void procs( const std::vector<EntityProc> & v ,
            const std::vector< Entity * > & m ,
            std::vector<unsigned> & p )
{
  typedef std::vector<EntityProc>::const_iterator Iter ;

  p.clear();

  std::vector<Iter> tmp ; tmp.reserve( m.size() );

  for ( unsigned j = 0 ; j < tmp.size() ; ++j ) {
    tmp.push_back( lower_bound( v , *m[j] ) );
  }

  if ( ! tmp.empty() ) {
    const Iter ve = v.end();
    Entity * const e0 = m[0] ;

    for ( ; tmp[0] != ve && tmp[0]->first == e0 ; ++tmp[0] ) {
      const unsigned proc = tmp[0]->second ;
      bool flag = true ;

      // Iterate remaining up to proc

      for ( unsigned j = 1 ; j < tmp.size() ; ++j ) {
        Entity * const ej = m[j] ;
        for ( ; tmp[j] != ve &&
                tmp[j]->first == ej &&
                tmp[j]->second < proc ; ++tmp[j] );
        if ( tmp[j] == ve ||
             tmp[j]->first != ej ||
             tmp[j]->second != proc ) {
          flag = false ;
        }
      }

      if ( flag ) { p.push_back( proc ); }
    }
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool verify( const std::vector<EntityProc> & v , std::string & msg )
{
  const LessEntityProc less_op ;

  bool result = true ;

  if ( ! v.empty() ) {

    const std::vector<EntityProc>::const_iterator e = v.end();
          std::vector<EntityProc>::const_iterator i ;

    std::ostringstream os ;

    for ( i = v.begin() ; result && i != e ; ++i ) {
      if ( ! ( result = i->first != NULL ) ) {
        msg.append( "Contains NULL entries" );
      }
    }

    if ( result ) {
      Mesh   & M = v[0].first->kernel().mesh();
      const Schema & S = M.schema();
      const unsigned p_size = S.parallel_size();

      for ( i = v.begin() ; result && i != e ; ++i ) {
        if ( ! ( result = & M == & i->first->kernel().mesh() ) ) {
          msg.append( "Contains entries from different meshes" );
        }
      }

      for ( i = v.begin() ; result && i != e ; ++i ) {
        const unsigned p = i->second ;
        if ( ! ( result = p < p_size ) ) {
          os << "Contains entry with bad processor " ;
          print_entry( os , *i );
          msg.append( os.str() );
        }
      }

      EntityProc old( v[0] ) ;

      for ( i = v.begin() ; result && ++i != e ; ) {
        if ( ! ( result = less_op( old , *i ) ) ) {
          os << "Contains out-of-order entries " ;
          print_entry( os , old );
          print_entry( os , *i );
          msg.append( os.str() );
        }
        old = *i ;
      }
    }
  }
  return result ;
}

bool comm_verify( ParallelMachine comm ,
                  const std::vector<EntityProc> & v ,
                  std::string & msg )
{
  static const char method[] = "phdmesh::comm_verify[symmetric]" ;

  // Verify local ordering, the result flag is parallel inconsistent

  bool result = verify( v , msg );

  // Verify parallel consistency
  // Communicate asymmetric and compare keys and owners

  const int p_rank = parallel_machine_rank( comm );
  const int p_size = parallel_machine_size( comm );

  std::ostringstream os ;

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( p_size , zero );

  const std::vector<EntityProc>::const_iterator i_end = v.end();
        std::vector<EntityProc>::const_iterator i ;

  if ( result ) {
    for ( i = v.begin() ; i != i_end ; ++i ) {
      send_size[ i->second ] += 2 * sizeof(unsigned long);
    }
  }

  CommAll comm_sparse ;
  {
    const unsigned * const p_send = & send_size[0] ;
    result = ! comm_sparse.allocate_buffers( comm , p_size / 4 ,
                                             p_send , NULL , ! result );
  }

  // The result flag is now parallel consistent

  if ( result ) {

    // Fill buffer
    for ( i = v.begin() ; i != i_end ; ++i ) {
      Entity & e = * i->first ;

      const unsigned long key   = e.key();
      const unsigned long owner = e.owner_rank();
      const unsigned      proc  = i->second ;

      CommBuffer & buf = comm_sparse.send_buffer( proc );
      buf.pack<unsigned long>( key );
      buf.pack<unsigned long>( owner );
    }

    comm_sparse.communicate();

    // Verify symmetry of sizes
    for ( int j = 0 ; result && j < p_size ; ++j ) {
      const unsigned nrecv = comm_sparse.recv_buffer( j ).remaining();
      if ( nrecv != send_size[j] ) {
        os << method ;
        os << " parallel inconsistency, P" << p_rank ;
        os << " expected from P" << j ;
        os << " " << ( send_size[j] / sizeof(unsigned long) ) ;
        os << " but received " ;
        os << ( nrecv / sizeof(unsigned long) ) ;
        os << " instead" ;
        result = false ;
        msg.append( os.str() );
      }
    }

    // Verify symmetry of content
    for ( i = v.begin() ; result && i != i_end ; ++i ) {
      Entity & e = * i->first ;

      const unsigned long this_key   = e.key();
      const unsigned long this_owner = e.owner_rank();
      const unsigned      proc  = i->second ;

      CommBuffer & buf = comm_sparse.recv_buffer( proc );
      unsigned long recv_key ;
      unsigned long recv_owner ;
      buf.unpack<unsigned long>( recv_key );
      buf.unpack<unsigned long>( recv_owner );

      if ( this_key != recv_key || this_owner != recv_owner ) {
        os << method ;
        os << " parallel inconsistency, P" << p_rank << " has " ;
        print_entity_key( os , this_key );
        os << "].owner(P" << this_owner ;
        os << ") versus " ;
        print_entity_key( os , recv_key );
        os << "].owner(P" << recv_owner ;
        os << ")" ;
        result = false ;
        msg.append( os.str() );
      }
    }

    // The result flag is now parallel inconsistent

    {
      unsigned flag = result ;
      all_reduce( comm , ReduceMin<1>( & flag ) );
      result = flag ;
    }
  }

  return result ;
}

//----------------------------------------------------------------------

bool comm_verify( ParallelMachine comm ,
                  const std::vector<EntityProc> & send ,
                  const std::vector<EntityProc> & recv ,
                  std::string & msg )
{
  static const char method[] = "phdmesh::comm_verify[asymmetric]" ;

  // Verify local ordering:

  bool result = verify( send , msg );

  const int p_rank = parallel_machine_rank( comm );
  const int p_size = parallel_machine_size( comm );

  std::ostringstream os ;

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( p_size , zero );
  std::vector<unsigned> recv_size( p_size , zero );

  std::vector<EntityProc>::const_iterator i ;

  if ( result ) {
    for ( i = send.begin() ; i != send.end() ; ++i ) {
      ++( send_size[ i->second ] );
    }
  }

  {
    unsigned msg_max = 0 ;
    unsigned * const p_send = & send_size[0] ;
    unsigned * const p_recv = & recv_size[0] ;
    result = ! comm_sizes( comm , p_size / 4 , msg_max ,
                           p_send , p_recv , ! result );
  }

  // Result flag is now parallel consistent

  if ( result ) {
    send_size.assign( p_size , zero );

    for ( i = recv.begin() ; i != recv.end() ; ++i ) {
      ++( send_size[ i->second ] );
    }

    for ( int j = 0 ; result && j < p_size ; ++j ) {
      const unsigned nrecv = recv_size[j] ;
      if ( nrecv != send_size[j] ) {
        os << method ;
        os << " parallel inconsistency, P" << p_rank ;
        os << " expected from P" << j ;
        os << " " << send_size[j] ;
        os << " but received " ;
        os << nrecv ;
        os << " instead" ;
        result = false ;
        msg.append( os.str() );
      }
    }

    // The result flag is now parallel inconsitent

    {
      unsigned flag = result ;
      all_reduce( comm , ReduceBitAnd<1>( & flag ) );
      result = flag ;
    }
  }
  return result ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_mesh_stats( Mesh & M ,
                      unsigned long * const counts ,
                      unsigned long * const max_id ,
                      bool local_flag )
{
  // Count locally owned entities

  const Schema & S = M.schema();
  ParallelMachine comm = S.parallel();
  Part & owns = S.owns_part();

  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    counts[i] = 0 ;
    max_id[i] = 0 ;

    const KernelSet & ks = M.kernels(  EntityType(i) );
    const EntitySet & es = M.entities( EntityType(i) );

    EntitySet::const_iterator ie ;
    KernelSet::const_iterator ik ;

    for ( ik = ks.begin() ; ik != ks.end() ; ++ik ) {
      if ( ik->has_superset( owns ) ) {
        counts[i] += ik->size();
      }
    }

    ie = es.end();
    if ( ie != es.begin() ) {
      --ie ;
      max_id[i] = ie->identifier();
    }
  }

  unsigned long flag = local_flag ? 1 : 0 ;

  all_reduce( comm ,
              ReduceSum<EntityTypeMaximum>( counts ) &
              ReduceMax<EntityTypeMaximum>( max_id ) &
              ReduceMax<1>( & flag ) );

  return flag ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

bool comm_mesh_entities(
  const EntityManager & manager ,
  Mesh & send_mesh ,
  Mesh & recv_mesh ,
  const std::vector<EntityProc> & send ,
        std::vector<EntityProc> & recv_info ,
  bool local_flag )
{
  enum { NUM_COUNT = EntityTypeMaximum + 1 };

  static const char method[] = "phdmesh::comm_mesh_entities" ;

  ParallelMachine comm = send_mesh.schema().parallel();

  if ( comm != recv_mesh.schema().parallel() ) {
    std::string msg ;
    msg.append( method );
    msg.append( " ERROR: given meshes with different parallel machines" );
    throw std::logic_error( msg );
  }

  CommAll all( comm );

  const unsigned u_zero = 0 ;
  const unsigned p_rank = all.parallel_rank();
  const unsigned p_size = all.parallel_size();

  const std::vector<EntityProc>::const_iterator i_end = send.end() ;
  const std::vector<EntityProc>::const_iterator i_beg = send.begin();
        std::vector<EntityProc>::const_iterator i ;

  // Verify monotonic ordering:
  {
    unsigned entity_type = 0 ;

    for ( i = i_beg ; i != i_end ; ++i ) {
      const unsigned etype = (unsigned) i->first->entity_type();

      if ( entity_type < etype ) { entity_type = etype ; }

      if ( etype != entity_type ) {
        std::ostringstream msg ;
        msg << "P" << p_rank << " " ;
        msg << method ;
        msg << " FAILED, send array is not monotonic" ;
        throw std::runtime_error( msg.str() );
      }
    }
  }

  std::vector<unsigned> count( p_size * NUM_COUNT , u_zero );

  // Sizing pass:
  {
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      all.send_buffer( p ).skip<unsigned>( NUM_COUNT );
    }

    for ( i = i_beg ; i != i_end ; ) {

      const std::vector<EntityProc>::const_iterator ib = i ;
      for ( ; i != send.end() && ib->first == i->first ; ++i );
      const std::vector<EntityProc>::const_iterator eb = i ;

      for ( i = ib ; i != eb ; ++i ) {
        if ( p_rank != i->second ) {
          CommBuffer & buf = all.send_buffer( i->second );

          const size_t size_before = buf.size();

          manager.send_entity( buf , recv_mesh , ib , eb );

          const size_t size_after = buf.size();

          if ( size_before < size_after ) {
            const unsigned etype  = (unsigned) i->first->entity_type();
            const unsigned offset = i->second * NUM_COUNT ;
            ++count[ offset ];
            ++count[ offset + etype + 1 ];
          }
        }
      }
    }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      if ( count[ p * NUM_COUNT ] == 0 ) {
        all.send_buffer( p ).reset();
      }
    }
  }

  // Allocate asymmetric send and receive buffers

  local_flag = all.allocate_buffers( p_size / 4 , local_flag );

  // Packing pass:
  {
    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      const unsigned * const tmp = & count[ p * NUM_COUNT ];
      if ( tmp[0] ) {
        CommBuffer & buf = all.send_buffer( p );
        buf.pack<unsigned>( tmp , NUM_COUNT );
      }
    }

    for ( i = i_beg ; i != i_end ; ) {

      const std::vector<EntityProc>::const_iterator ib = i ;
      for ( ; i != send.end() && ib->first == i->first ; ++i );
      const std::vector<EntityProc>::const_iterator eb = i ;

      for ( i = ib ; i != eb ; ++i ) {
        if ( p_rank != i->second ) {
          CommBuffer & buf = all.send_buffer( i->second );
          manager.send_entity( buf , recv_mesh , ib , eb );
        }
      }
    }
  }

  all.communicate();

  // Unpacking pass:
  {
    for ( unsigned recv_p = 0 ; recv_p < p_size ; ++recv_p ) {
      CommBuffer & buf = all.recv_buffer( recv_p );
      unsigned * const tmp = & count[ recv_p * NUM_COUNT ];
      Copy<NUM_COUNT>( tmp , u_zero );
      if ( buf.capacity() ) {
        buf.unpack<unsigned>( tmp , NUM_COUNT );
      }
    }

    for ( unsigned etype = 0 ; etype < EntityTypeMaximum ; ++etype ) {
      for ( unsigned recv_p = 0 ; recv_p < p_size ; ++recv_p ) {
        const unsigned num = count[ recv_p * NUM_COUNT + etype + 1 ];
        CommBuffer & buf = all.recv_buffer( recv_p );
        for ( unsigned j = 0 ; j < num ; ++j ) {
          manager.receive_entity( buf , recv_mesh , recv_p , recv_info );
        }
      }
    }
  }

  return local_flag ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Heterogeneity?

bool comm_mesh_field_values(
  const std::vector<const Field<void,0> *> & fields ,
  const std::vector<EntityProc> & domain ,
  const std::vector<EntityProc> & range ,
  bool local_flag )
{
  if ( fields.empty() ) { return local_flag ; }

  const Schema & schema = fields[0]->schema();

  const unsigned parallel_size = schema.parallel_size();
  const unsigned parallel_rank = schema.parallel_rank();
  const bool     asymmetric    = & domain != & range ;

  // Memory required per field

  unsigned fields_size[ EntityTypeMaximum ] ;
  for ( unsigned i = 0 ; i < EntityTypeMaximum ; ++i ) {
    fields_size[i] = 0 ;
  }

  const std::vector<const Field<void,0> *>::const_iterator fe = fields.end();
  const std::vector<const Field<void,0> *>::const_iterator fb = fields.begin();
        std::vector<const Field<void,0> *>::const_iterator fi ;

  for ( fi = fb ; fi != fe ; ++fi ) {
    const Field<void,0> & f   = **fi ;
    const EntityType      et  = f.entity_type();
    const unsigned        fsize = field_type_size( f.field_types_ordinal() );

    fields_size[ et ] += f.max_length() * fsize ;
  }

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  std::vector<EntityProc>::const_iterator i ;

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const Field<void,0> & f = **fi ;
        const EntityType      fet = f.entity_type();
        const EntityType      eet = e.entity_type();
        if ( fet == eet ) {
          const FieldDimension & dim = dimension( f , e );
          e_size += dim.size();
        }
      } 
      send_size[ p ] += e_size ;
    }
  }

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      unsigned e_size = 0 ;
      for ( fi = fb ; fi != fe ; ++fi ) {
        const Field<void,0> & f = **fi ;
        const EntityType      fet = f.entity_type();
        const EntityType      eet = e.entity_type();
        if ( fet == eet ) {
          const FieldDimension & dim = dimension( f , e );
          e_size += dim.size();
        }
      } 
      recv_size[ p ] += e_size ;
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const s_size = & send_size[0] ;
    const unsigned * const r_size = & recv_size[0] ;
    sparse.allocate_buffers( schema.parallel(), parallel_size / 4 , s_size, r_size);
  }

  // Pack for send:

  for ( i = domain.begin() ; i != domain.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || parallel_rank == e.owner_rank() ) {
      CommBuffer & b = sparse.send_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const Field<void,0> & f = **fi ;
        const EntityType      fet = f.entity_type();
        const EntityType      eet = e.entity_type();
        if ( fet == eet ) {
          const FieldDimension & dim = dimension( f , e );
          unsigned char * data = (unsigned char *) e.data( f );
          b.pack<unsigned char>( data , dim.size() );
        }
      }
    }
  }

  // Communicate:

  sparse.communicate();

  // Unpack for recv:

  for ( i = range.begin() ; i != range.end() ; ++i ) {
    Entity       & e = * i->first ;
    const unsigned p = i->second ;

    if ( asymmetric || p == e.owner_rank() ) {
      CommBuffer & b = sparse.recv_buffer( p );
      for ( fi = fb ; fi != fe ; ++fi ) {
        const Field<void,0> & f = **fi ;
        const EntityType      fet = f.entity_type();
        const EntityType      eet = e.entity_type();
        if ( fet == eet ) {
          const FieldDimension & dim = dimension( f , e );
          unsigned char * data = (unsigned char *) e.data( f );
          b.unpack<unsigned char>( data , dim.size() );
        }
      }
    }
  }

  return local_flag ;
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

}

