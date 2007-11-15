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

#include <math.h>
#include <sstream>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <util/OctTreeOps.hpp>

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>
#include <mesh/Proximity.hpp>

#include <mesh_io/ExoII.hpp>

#include "Gears.hpp"

using namespace phdmesh ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

#if defined( PHDMESH_HAS_MPI )

void parallel_gather( ParallelMachine comm ,
                      unsigned        rank ,
                      const unsigned * local ,
                      unsigned * global ,
                      unsigned count )
{
  unsigned * tmp = const_cast<unsigned*>( local );

  MPI_Gather( tmp , count , MPI_UNSIGNED ,
              global , count , MPI_UNSIGNED , rank , comm );
}

#else

void parallel_gather( ParallelMachine ,
                      unsigned ,
                      const unsigned * local ,
                      unsigned * global ,
                      unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

#endif

}

//----------------------------------------------------------------------

void test_six_gears_face_proximity(
  Mesh & M ,
  const Field<double,1> & gear_coordinates ,
  const Field<double,1> & field_proximity ,
  const ProximitySearch & prox_search ,
  EntityProcSet & domain ,
  EntityProcSet & range )
{
  static const char method[] = "phdmesh::test_six_gears_face_proximity" ;

  const Schema & schema = M.schema();
  const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS failed parallel consistency before proximity"
                << std::endl ;
    }
    throw std::runtime_error(std::string("SIX_GEARS"));
  }

  domain.clear();
  range.clear();

  std::vector< std::pair<IdentProc,IdentProc> > proximity ;

  double dt_proximity_search = wall_time();

  const unsigned num_search_tasks =
    proximity_search( M , prox_search , Face , proximity );

  dt_proximity_search = wall_time() - dt_proximity_search ;

  const std::vector< std::pair<IdentProc,IdentProc> >::iterator
      i_end = proximity.end() ,
      i_beg = proximity.begin() ;
  std::vector< std::pair<IdentProc,IdentProc> >::iterator i ;

  // Set "in proximity" flags on nodes in a contact face
  {
    // Clear the node value
    const KernelSet::const_iterator k_beg = M.kernels( Node ).begin();
    const KernelSet::const_iterator k_end = M.kernels( Node ).end();
    KernelSet::const_iterator k ;
    for ( k = k_beg ; k != k_end ; ++k ) {
      unsigned n = k->size();
      double * const data = k->data( field_proximity );
      double * const coord = k->data( gear_coordinates );
      for ( unsigned j = 0 ; j < n ; ++j ) { data[j] = coord[1+j*3] ; }
    }

    for ( i = i_beg ; i != i_end ; ++i ) {
      IdentProc d[2] ;
      d[0] = i->first ;
      d[1] = i->second ;

      for ( unsigned j = 0 ; j < 2 ; ++j ) {
        if ( p_rank == d[j].proc ) {
          Entity & face = * M.get_entity( Face , d[j].ident , method );
          for ( ConnectSpan face_nodes = face.connections( Node , Uses );
                face_nodes ; ++face_nodes ) {
            Entity & node = * face_nodes->entity();
            double * const data = node.data( field_proximity );
            *data = 20 ;
          }
        }
      }
    }
  }

  // Count number of contact faces
  {
    unsigned d_count = 0 ;
    unsigned r_count = 0 ;

    for ( i = i_beg ; i_end != i ; ++i ) {
      const IdentProc d = i->first ;
      const IdentProc r = i->second ;

      if ( d.proc == p_rank ) { ++d_count ; }
      if ( r.proc == p_rank ) { ++r_count ; }
    }

    const unsigned u_zero = 0 ;
    std::vector<unsigned> recv_counts( p_size * 2 , u_zero );

    Part & owns_part = schema.owns_part();

    unsigned send_count[2] = { 0 , 0 };
    unsigned long counts[ EntityTypeMaximum ];

    partset_entity_count( M , owns_part , counts );

    send_count[0] = counts[2] ; // Owned faces
    send_count[1] = d_count ;   // Owned domain proximity face

    parallel_gather( M.parallel() , 0 ,
                     send_count , & recv_counts[0] , 2 );

    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS surface proximity search results" << std::endl ;
      unsigned total = 0 ;
      for ( unsigned p = 0 ; p < p_size ; ++p ) {
        std::cout << "  P" << p << " has " << recv_counts[2*p+1] 
                  << " domain-owned relations" << std::endl ;
        total += recv_counts[2*p+1];
      }
      std::cout << "  Total of " << total << " domain-owned relations"
                << std::endl ;
      std::cout << "  Proximity search used "
                << dt_proximity_search
                << " seconds with "
                << num_search_tasks
                << " parallel tasks"
                << std::endl ;
      std::cout.flush();
    }
  }

  // Create an 'Other' entity to connect the domain entity to the range entity
  // This entity will be owned by the domain processor will appear in the
  // 'aura' of the range processor.
  // If the range face is not on the domain processor then share it.

  EntityProcSet to_be_shared ;

  for ( i = i_beg ; i_end != i ; ++i ) {
    const IdentProc & d = i->first ;
    const IdentProc & r = i->second ;

    if ( r.proc == p_rank && d.proc != p_rank ) {

      // The range face needs to be shared with the domain processor

      EntityProc ep ;
      ep.first  = M.get_entity( Face , r.ident , method );
      ep.second = d.proc ;
      to_be_shared.push_back( ep );

      for ( ConnectSpan con = ep.first->connections(); con ; ++con ) {
        if ( con->type() == Uses ) {
          ep.first = con->entity();
          to_be_shared.push_back( ep );
        }
      }
    }
  }

  sort_unique( to_be_shared );

  const EntityProcSet sharing_A( M.shares() );

  comm_mesh_add_sharing( M , to_be_shared );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS failed parallel consistency of add_sharing"
                << std::endl ;
    }
    throw std::runtime_error(std::string("SIX_GEARS"));
  }

  comm_mesh_scrub_sharing( M );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS failed parallel consistency of scrub_sharing"
                << std::endl ;
    }
    throw std::runtime_error(std::string("SIX_GEARS"));
  }

  unsigned flag = M.shares() == sharing_A ;

  all_reduce( M.parallel() , ReduceMin<1>( & flag ) );

  if ( ! flag ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS failed parallel sharing before == after"
                << std::endl ;
    }
    throw std::runtime_error(std::string("SIX_GEARS"));
  }
}

//----------------------------------------------------------------------

void test_six_gears( ParallelMachine pm , std::istream & )
{
  const double TWO_PI = 2.0 * acos( (double) -1.0 );

  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );

  const unsigned kernel_capacity[ EntityTypeMaximum ] =
    { 100 , 100 , 100 , 100 , 100 };
    // { 20 , 20 , 20 , 20 , 20 };

  Schema S ;

  GearFields gear_fields( S );

  //------------------------------
  // Six circular gears with a 'z' axis of rotation.
  // They are spaced every 60 degrees about the origin

  const double sqrt_3 = sqrt( 3.0 );
  const double center_A[3] = {  2.0 ,      0.0 , 0.0 };
  const double center_B[3] = {  1.0 ,   sqrt_3 , 0.0 };
  const double center_C[3] = { -1.0 ,   sqrt_3 , 0.0 };
  const double center_D[3] = { -2.0 ,      0.0 , 0.0 };
  const double center_E[3] = { -1.0 , - sqrt_3 , 0.0 };
  const double center_F[3] = {  1.0 , - sqrt_3 , 0.0 };

  // Exactly touch = 1.0, force overlap by adding a little
  const double rad_max = 1.0 + 0.001 ;
  const double rad_min = 0.6 ;
  const double z_min   = -0.4 ;
  const double z_max   =  0.4 ;

  const double elem_h = 0.10 ;

  const unsigned angle_num = (unsigned) ( TWO_PI / elem_h );
  const unsigned rad_num   = (unsigned) ( 1 + ( rad_max - rad_min ) / elem_h );
  const unsigned z_num     = (unsigned) ( 1 + ( z_max   - z_min )   / elem_h );
  const unsigned tot_num   = angle_num * ( rad_num - 1 ) * ( z_num - 1 );

  if ( p_rank == 0 ) {
    std::cout << std::endl
              << "SIX_GEARS meshing: #Angles = " << angle_num
              << " , #Radial = " << rad_num
              << " , #Thickness = " << z_num
              << " , Elements/gear = " << tot_num
              << std::endl ;
  }
  
  Gear A( S , std::string("A") , gear_fields ,
          center_A ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , 1 );

  Gear B( S , std::string("B") , gear_fields ,
          center_B ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , -1 );

  Gear C( S , std::string("C") , gear_fields ,
          center_C ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , 1 );

  Gear D( S , std::string("D") , gear_fields ,
          center_D ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , -1 );

  Gear E( S , std::string("E") , gear_fields ,
          center_E ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , 1 );

  Gear F( S , std::string("F") , gear_fields ,
          center_F ,
          rad_min , rad_max , rad_num ,
          z_min , z_max , z_num ,
          angle_num , -1 );

#if 0
  std::vector<Gear*> gears( i_end * j_end );

  for ( unsigned j = 0 ; j < j_end ; ++j ) {
    double y = j * sqrt_3 ;
    double x ;
    for ( unsigned i = 0 ; i < i_end ; ++i ) {
      int dir = i % 2 ? 1 : -1 ;

      if ( j % 2 ) { // Odd
        x = i * 3 + i % 2 ;
      }
      else { // Even
        x = i * 3 + ( 1 - i % 2 );
      }

      std::ostream name ; name << "G_" << i << "_" << j ;

      Gear * g = new Gear( S , name.str() , gear_fields ,
                           rad_min , rad_max , rad_num ,
                           z_min , z_max , z_num ,
                           angle_num , dir );

      gears[ j * i_end + i ] = g ;

      S.declare_part_attribute( g->m_surf , & proximity_search , NULL );

      file_schema.declare_part( g->m_gear.name() , 1 , exodus::HEX , 8 );
    }
  }

  S.commit();
#endif

  //------------------------------
  // Proximity search.
  // Tag the surface parts with the proximity search object.
  // This prevents self-search of faces within a single gear.

  ProximitySearch proximity_search( gear_fields.current_coord , 0.25 );

  S.declare_part_attribute( A.m_surf , & proximity_search , NULL );
  S.declare_part_attribute( B.m_surf , & proximity_search , NULL );
  S.declare_part_attribute( C.m_surf , & proximity_search , NULL );
  S.declare_part_attribute( D.m_surf , & proximity_search , NULL );
  S.declare_part_attribute( E.m_surf , & proximity_search , NULL );
  S.declare_part_attribute( F.m_surf , & proximity_search , NULL );

  //------------------------------

  exodus::FileSchema file_schema( S , gear_fields.model_coord , 
                                      gear_fields.element_attr );

  file_schema.declare_part( A.m_gear.name() , 1 , exodus::HEX , 8 );
  file_schema.declare_part( B.m_gear.name() , 2 , exodus::HEX , 8 );
  file_schema.declare_part( C.m_gear.name() , 3 , exodus::HEX , 8 );
  file_schema.declare_part( D.m_gear.name() , 4 , exodus::HEX , 8 );
  file_schema.declare_part( E.m_gear.name() , 5 , exodus::HEX , 8 );
  file_schema.declare_part( F.m_gear.name() , 6 , exodus::HEX , 8 );

  Field<double,1> & field_node_proximity =
    S.declare_field<double,1>( Node , std::string("proximity") , 1 );

  S.declare_field_dimension( field_node_proximity , S.universal_part() , 1 );

  S.commit();

  Mesh M(S,pm,kernel_capacity);

  A.mesh(M);
  B.mesh(M);
  C.mesh(M);
  D.mesh(M);
  E.mesh(M);
  F.mesh(M);

  comm_mesh_discover_sharing( M );
  comm_mesh_regenerate_aura( M );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    std::cout << "SIX_GEARS Failed parallel consistency" << std::endl ;
    return ;
  }

  // Copy coordinates to the aura nodes
  {
    std::vector< const Field<void,0> *> fields ;
    const Field<void,0> * ptr = NULL ;

    ptr = & gear_fields.gear_coord ;    fields.push_back( ptr );
    ptr = & gear_fields.model_coord ;   fields.push_back( ptr );
    ptr = & gear_fields.current_coord ; fields.push_back( ptr );

    comm_mesh_field_values( M, M.aura_domain(), M.aura_range(), fields, false );
  }

  {
    unsigned long counts[ EntityTypeMaximum ];
    unsigned long max_id[ EntityTypeMaximum ];

    comm_mesh_stats( M , counts , max_id );

    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS Global Counts = " 
                << "node = " << counts[0] << " "
                << "edge = " << counts[1] << " "
                << "face = " << counts[2] << " "
                << "elem = " << counts[3] << " "
                << "other = " << counts[4] << std::endl ;
      std::cout << "SIX_GEARS Global MaxId  = " 
                << "node = " << max_id[0] << " "
                << "edge = " << max_id[1] << " "
                << "face = " << max_id[2] << " "
                << "elem = " << max_id[3] << " "
                << "other = " << max_id[4] << std::endl ;
    }
  }


  if ( ! comm_verify_shared_entity_values( M , gear_fields.gear_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.gear_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  if ( ! comm_verify_shared_entity_values( M , gear_fields.model_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "TWO_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.model_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  if ( ! comm_verify_shared_entity_values( M , gear_fields.current_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "SIX_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.current_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  exodus::FileOutput * exo = NULL ;

  {
    file_schema.assign_indices( M );

    std::ostringstream file_name ;
    file_name << "six_gears_np" << p_size << ".exo" ;

    std::string title( "PHDMESH Gears test problem" );

    std::vector< const Field<void,0> * > out_fields ;
    const Field<void,0> * tmp ;

    // tmp = & gear_fields.gear_coord ;    out_fields.push_back( tmp );
    // tmp = & gear_fields.current_coord ; out_fields.push_back( tmp );
    tmp = & gear_fields.displacement ;  out_fields.push_back( tmp );
    tmp = & field_node_proximity ;      out_fields.push_back( tmp );

    int flags[ EntityTypeMaximum ] = { 0 , 0 , 0 , 1 , 0 };

    double dt = wall_time();
    exo = new exodus::FileOutput( file_schema, M,
                                  file_name.str(), title,
                                  false , out_fields, flags);
    dt = wall_time() - dt ;
    if ( p_rank == 0 ) {
      std::cout << "Exodus model output = " << dt << " seconds" << std::endl ;
    }
  }

  EntityProcSet prox_domain ;
  EntityProcSet prox_range ;

  test_six_gears_face_proximity( M ,
                                 gear_fields.gear_coord ,
                                 field_node_proximity ,
                                 proximity_search ,
                                 prox_domain , prox_range );

  std::vector<OctTreeKey> rebal_cut_keys ;

  comm_mesh_rebalance( M , gear_fields.current_coord , NULL , rebal_cut_keys );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    std::cout << "SIX_GEARS Failed parallel rebalance consistency"
              << std::endl ;
    return ;
  }

  if ( p_rank == 0 ) {
    std::cout << "SIX_GEARS Rebalance successful" << std::endl ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      std::cout << "  P" << i << ": cut = " << rebal_cut_keys[i] << std::endl ;
    }

    std::cout.flush();
  }

  //------------------------------
  // Turn the gears in opposite directions by the same amount

  {
    const unsigned nsteps = 120 ;

    for ( unsigned i = 0 ; i <= nsteps ; ++i ) {

      if ( p_rank == 0 ) {
        std::cout << std::endl << "SIX_GEARS turn "
                  << i << " of " << nsteps << std::endl ;
      }

      M.update_state();

      const double angle = ( TWO_PI * i ) / ( (double) nsteps );

      // Set coordinates[ STATE_NEW ] to turned coordinates

      A.turn( angle );
      B.turn( angle );
      C.turn( angle );
      D.turn( angle );
      E.turn( angle );
      F.turn( angle );

      // Copy the coordinates to the aura nodes
      {
        std::vector< const Field<void,0> *> fields ;
        const Field<void,0> * const ptr = & gear_fields.current_coord ;
        fields.push_back( ptr );
        comm_mesh_field_values( M, M.aura_domain(),
                                   M.aura_range(), fields, false );
      }

      // Check parallel consistency of shared variable

      if ( ! comm_verify_shared_entity_values(M,gear_fields.current_coord) ) {
        if ( p_rank == 0 ) {
          std::cout << "SIX_GEARS FAILED for shared values of " ;
          std::cout << gear_fields.current_coord.name();
          std::cout << std::endl ;
        }
        return ;
      }

      const Field<double,1> & node_current_coord_old =
        gear_fields.current_coord[ StateOld ];

      if ( ! comm_verify_shared_entity_values(M , node_current_coord_old) ) {
        if ( p_rank == 0 ) {
          std::cout << "SIX_GEARS FAILED for shared values of " ;
          std::cout << node_current_coord_old.name();
          std::cout << std::endl ;
        }
        return ;
      }

      test_six_gears_face_proximity( M ,
                                     gear_fields.gear_coord ,
                                     field_node_proximity ,
                                     proximity_search ,
                                     prox_domain , prox_range );

      // Surface in proximity are added to the Aura




      {
        double dt = wall_time();
        exo->write( 0.0 );
        dt = wall_time() - dt ;
        if ( p_rank == 0 ) {
          std::cout << "Exodus model results = "
                    << dt << " seconds" << std::endl ;
        }
      }
    }
  }

  delete exo ; exo = NULL ;

  //------------------------------
  // Test reading the written file
  {
    Schema S_read ;
    
    GearFields gear_fields_read( S_read );

    std::ostringstream i_file_name , o_file_name ;
    i_file_name << "six_gears_np" << p_size << ".exo" ;
    o_file_name << "six_gears_io_np" << p_size << ".exo" ;

    exodus::FileSchema FS_read( S_read , gear_fields_read.model_coord , 
                                         gear_fields_read.element_attr ,
                                         i_file_name.str() ,
                                         pm );

    S_read.commit();

    Mesh M_read(S_read,pm,kernel_capacity);

    std::vector< const Field<void,0> * > none ;

    exodus::FileInput FI_read( FS_read , M_read , i_file_name.str() , none );

    std::string title( "TEST INPUT/OUTPUT" );

    exodus::FileOutput FO_read(FS_read,M_read,
                                o_file_name.str(),title,true,none);
    FO_read.write(0.0);
  }

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << "SIX_GEARS test successful" << std::endl ;
    std::cout.flush();
  }
}

