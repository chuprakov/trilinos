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

#include <util/TPI.h>
#include <util/ParallelComm.hpp>
#include <util/ParallelReduce.hpp>
#include <util/OctTreeOps.hpp>

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/FieldData.hpp>
#include <mesh/Comm.hpp>
#include <mesh/Proximity.hpp>

#include <mesh_io/ExoII.hpp>

#include "Gears.hpp"

using namespace phdmesh ;

typedef GearFields::CylindricalField CylindricalField ;
typedef GearFields::CartesianField   CartesianField ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

#if 0

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

#endif

}

//----------------------------------------------------------------------

void test_gears_face_proximity(
  Mesh & M ,
  const CylindricalField & gear_coordinates ,
  const Field<double>  & field_proximity ,
  const ProximitySearch & prox_search ,
  EntityProcSet & domain ,
  EntityProcSet & range ,
  const bool verify ,
  double & dt_proximity_search ,
  double & dt_ghosting )
{
  static const char method[] = "phdmesh::test_gears_face_proximity" ;

  // const unsigned p_size = M.parallel_size();
  const unsigned p_rank = M.parallel_rank();
  double wt ;

  if ( verify && ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS failed parallel consistency before proximity"
                << std::endl ;
    }
    throw std::runtime_error(std::string("N_GEARS"));
  }

  domain.clear();
  range.clear();

  std::vector< std::pair<IdentProc,IdentProc> > proximity ;

  wt = wall_time();

  proximity_search( M , prox_search , Face , proximity );

  dt_proximity_search += wall_dtime(wt);

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
      double * const data  = field_data( field_proximity , *k );
      double * const coord = field_data( gear_coordinates , *k );
      for ( unsigned j = 0 ; j < n ; ++j ) { data[j] = coord[1+j*3] ; }
    }

    for ( i = i_beg ; i != i_end ; ++i ) {
      IdentProc d[2] ;
      d[0] = i->first ;
      d[1] = i->second ;

      for ( unsigned j = 0 ; j < 2 ; ++j ) {
        if ( p_rank == d[j].proc ) {
          const entity_key_type key = entity_key( Face , d[j].ident );
          Entity & face = * M.get_entity( key , method );
          for ( RelationSpan face_nodes = face.relations( Node );
                face_nodes ; ++face_nodes ) {
            Entity & node = * face_nodes->entity();
            double * const data = field_data( field_proximity , node );
            *data = 20 ;
          }
        }
      }
    }
  }

  dt_ghosting += wall_dtime(wt);

  wt = wall_time();

  // Create an 'Other' entity to relation the domain entity to the range entity
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
      ep.first  = M.get_entity( entity_key( Face , r.ident ) , method );
      ep.second = d.proc ;
      to_be_shared.push_back( ep );

      for ( RelationSpan con = ep.first->relations(); con ; ++con ) {
        if ( con->entity_type() < Face ) {
          ep.first = con->entity();
          to_be_shared.push_back( ep );
        }
      }
    }
  }

  sort_unique( to_be_shared );

  const EntityProcSet sharing_A( M.shares() );

  comm_mesh_add_sharing( M , to_be_shared );

  dt_ghosting += wall_dtime( wt );

  if ( verify && ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS failed parallel consistency of add_sharing"
                << std::endl ;
    }
    throw std::runtime_error(std::string("N_GEARS"));
  }

  comm_mesh_scrub_sharing( M );

  if ( verify && ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS failed parallel consistency of scrub_sharing"
                << std::endl ;
    }
    throw std::runtime_error(std::string("N_GEARS"));
  }

  unsigned flag = M.shares() == sharing_A ;

  all_reduce( M.parallel() , ReduceMin<1>( & flag ) );

  if ( ! flag ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS failed parallel sharing before == after"
                << std::endl ;
    }
    throw std::runtime_error(std::string("N_GEARS"));
  }
}

//----------------------------------------------------------------------

void test_gears( ParallelMachine pm ,
                 const unsigned i_end ,
                 const unsigned j_end ,
                 const unsigned k_end ,
                 const std::string & exo_file_name ,
                 const bool verify )
{
  const double TWO_PI = 2.0 * acos( (double) -1.0 );

  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );

  const unsigned kernel_capacity = 100 ; // 20 ;

  Schema S ;

  GearFields gear_fields( S );

  double dt_rebalance = 0 ;
  double dt_proximity = 0 ;
  double dt_ghosting = 0 ;
  double dt_exo_write = 0 ;

  //------------------------------
  // Six circular gears with a 'z' axis of rotation.
  // They are spaced every 60 degrees about the origin

  const double sqrt_3 = sqrt( 3.0 );

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
              << "N_GEARS meshing: #Angles = " << angle_num
              << " , #Radial = " << rad_num
              << " , #Thickness = " << z_num
              << " , Elements/gear = " << tot_num
              << std::endl ;
  }
  
  //------------------------------
  // Proximity search.
  // Tag the surface parts with the proximity search object.
  // This prevents self-search of faces within a single gear.

  ProximitySearch proximity_search( gear_fields.current_coord , 0.25 );

  Field<double> & field_node_proximity =
    S.declare_field< Field<double> >( std::string("proximity") );

  S.declare_field_exists( field_node_proximity , Node , S.universal_part() );

  //------------------------------

  std::vector<Gear*> gears( i_end * j_end * k_end );

  for ( unsigned k = 0 ; k < k_end ; ++k ) {
    for ( unsigned j = 0 ; j < j_end ; ++j ) {
      double center[3] ;
      center[2] = k - z_min ;
      center[1] = sqrt_3 * j ;
      for ( unsigned i = 0 ; i < i_end ; ++i ) {
        int dir = i % 2 ? 1 : -1 ;

        if ( j % 2 ) { // Odd
          center[0] = i * 3 + i % 2 ;
          dir = - dir ;
        }
        else { // Even
          center[0] = i * 3 + ( 1 - i % 2 );
        }


        std::ostringstream name ; name << "G_" << i << "_" << j << "_" << k ;

        Gear * g = new Gear( S , name.str() , gear_fields ,
                             center ,
                             rad_min , rad_max , rad_num ,
                             z_min , z_max , z_num ,
                             angle_num , dir );

        gears[ k * j_end * i_end + j * i_end + i ] = g ;

        S.declare_part_attribute<ProximitySearch>( g->m_surf ,
                                                   & proximity_search ,
                                                   false );
      }
    }
  }

  //------------------------------

  exodus::FileSchema file_schema( S , gear_fields.model_coord , 
                                      gear_fields.element_attr );

  {
    unsigned j = 1 ;
    for ( std::vector<Gear*>::iterator
          i = gears.begin() ; i != gears.end() ; ++i , ++j ) {
      file_schema.declare_part( (*i)->m_gear , j , exodus::HEX , 8 );
    }
  }

  S.commit();

  //------------------------------

  Mesh M(S,pm,kernel_capacity);

  double wt = wall_time();

  for ( std::vector<Gear*>::iterator
        i = gears.begin() ; i != gears.end() ; ++i ) {
    (*i)->mesh( M );
  }

  comm_mesh_discover_sharing( M );
  comm_mesh_regenerate_aura( M );

  double dt_mesh_gen = wall_dtime( wt );

  if ( verify && ! comm_mesh_verify_parallel_consistency( M ) ) {
    std::cout << "N_GEARS Failed parallel consistency" << std::endl ;
    return ;
  }

  // Copy coordinates to the aura nodes
  {
    std::vector< const FieldBase *> fields ;
    const FieldBase * ptr = NULL ;

    ptr = & gear_fields.gear_coord ;    fields.push_back( ptr );
    ptr = & gear_fields.model_coord ;   fields.push_back( ptr );
    ptr = & gear_fields.current_coord ; fields.push_back( ptr );

    comm_mesh_field_values( M, M.aura_domain(), M.aura_range(), fields, false );
  }

  {
    entity_id_type counts[ EntityTypeEnd ];
    entity_id_type max_id[ EntityTypeEnd ];

    comm_mesh_stats( M , counts , max_id );

    if ( p_rank == 0 ) {
      std::cout << "N_GEARS Global Counts { " 
                << "node = " << counts[0] << " , "
                << "edge = " << counts[1] << " , "
                << "face = " << counts[2] << " , "
                << "elem = " << counts[3] << " , "
                << "particle = " << counts[4] << " , "
                << "constraint = " << counts[5] << " }" << std::endl ;
      std::cout << "N_GEARS Global MaxId { " 
                << "node = " << max_id[0] << " , "
                << "edge = " << max_id[1] << " , "
                << "face = " << max_id[2] << " , "
                << "elem = " << max_id[3] << " , "
                << "particle = " << max_id[4] << " , "
                << "constraint = " << max_id[5] << " }" << std::endl ;
    }
  }


  if ( verify && ! comm_verify_shared_entity_values( M , Node , gear_fields.gear_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.gear_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  if ( verify && ! comm_verify_shared_entity_values( M , Node , gear_fields.model_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "TWO_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.model_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  if ( verify && ! comm_verify_shared_entity_values( M , Node , gear_fields.current_coord ) ) {
    if ( p_rank == 0 ) {
      std::cout << "N_GEARS FAILED for shared values of " ;
      std::cout << gear_fields.current_coord.name();
      std::cout << std::endl ;
    }
    return ;
  }

  exodus::FileOutput * exo = NULL ;

  if ( exo_file_name.size() ) {
    file_schema.assign_indices( M );

    std::string title( "PHDMESH Gears test problem" );

    std::vector< const FieldBase * > out_fields ;
    const FieldBase * tmp ;

    // tmp = & gear_fields.gear_coord ;    out_fields.push_back( tmp );
    // tmp = & gear_fields.current_coord ; out_fields.push_back( tmp );
    tmp = & gear_fields.displacement ;  out_fields.push_back( tmp );
    tmp = & field_node_proximity ;      out_fields.push_back( tmp );

    int flags[ EntityTypeEnd ] = { 0 , 0 , 0 , 1 , 0 , 0 };

    wt = wall_time();
    exo = new exodus::FileOutput( file_schema, M,
                                  exo_file_name , title,
                                  false , out_fields, flags);
    dt_exo_write += wall_dtime(wt);
  }

  EntityProcSet prox_domain ;
  EntityProcSet prox_range ;

  test_gears_face_proximity( M ,
                             gear_fields.gear_coord ,
                             field_node_proximity ,
                             proximity_search ,
                             prox_domain , prox_range ,
                             verify ,
                             dt_proximity , dt_ghosting );

  wt = wall_time();

  std::vector<OctTreeKey> rebal_cut_keys ;

  comm_mesh_rebalance( M , gear_fields.current_coord , NULL , rebal_cut_keys );

  dt_rebalance = wall_dtime( wt );

  if ( verify && ! comm_mesh_verify_parallel_consistency( M ) ) {
    std::cout << "N_GEARS Failed parallel rebalance consistency"
              << std::endl ;
    return ;
  }

  if ( p_rank == 0 ) {
    std::cout << "N_GEARS Rebalance successful" << std::endl ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      std::cout << "  P" << i << ": cut = " << rebal_cut_keys[i] << std::endl ;
    }

    std::cout.flush();
  }

  //------------------------------
  // Turn the gears in opposite directions by the same amount

  dt_proximity = 0 ;
  dt_ghosting = 0 ;

  const unsigned nsteps = 120 ;

  {
    for ( unsigned i = 0 ; i <= nsteps ; ++i ) {

      M.update_state();

      const double angle = ( TWO_PI * i ) / ( (double) nsteps );

      // Set coordinates[ STATE_NEW ] to turned coordinates

      for ( std::vector<Gear*>::iterator
            j = gears.begin() ; j != gears.end() ; ++j ) {
        (*j)->turn( angle );
      }

      // Copy the coordinates to the aura nodes
      {
        std::vector< const FieldBase *> fields ;
        const FieldBase * const ptr = & gear_fields.current_coord ;
        fields.push_back( ptr );
        comm_mesh_field_values( M, M.aura_domain(),
                                   M.aura_range(), fields, false );
      }

      // Check parallel consistency of shared variable

      if ( verify && ! comm_verify_shared_entity_values(M,Node,gear_fields.current_coord) ) {
        if ( p_rank == 0 ) {
          std::cout << "N_GEARS FAILED for shared values of " ;
          std::cout << gear_fields.current_coord.name();
          std::cout << std::endl ;
        }
        return ;
      }

      const CartesianField & node_current_coord_old =
        gear_fields.current_coord[ StateOld ];

      if ( verify && ! comm_verify_shared_entity_values(M , Node, node_current_coord_old) ) {
        if ( p_rank == 0 ) {
          std::cout << "N_GEARS FAILED for shared values of " ;
          std::cout << node_current_coord_old.name();
          std::cout << std::endl ;
        }
        return ;
      }

      test_gears_face_proximity( M ,
                                 gear_fields.gear_coord ,
                                 field_node_proximity ,
                                 proximity_search ,
                                 prox_domain , prox_range ,
                                 verify ,
                                 dt_proximity , dt_ghosting );

      // Surface in proximity are added to the Aura




      if ( NULL != exo ) {
        wt = wall_time();
        exo->write( 0.0 );
        dt_exo_write += wall_dtime( wt );
      }
    }
  }

  if ( NULL != exo ) { delete exo ; exo = NULL ; }

  {
    for ( std::vector<Gear*>::iterator
          i = gears.begin() ; i != gears.end() ; ++i ) {
      delete *i ;
    }
    gears.clear();
  }

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << "N_GEARS " << nsteps
              << ", mesh = " << dt_mesh_gen 
              << ", rebal = " << dt_rebalance 
              << ", search = " << dt_proximity 
              << ", ghost = " << dt_ghosting
              << ", write = " << dt_exo_write
              << std::endl ;
    std::cout.flush();
  }
}

void test_gears( phdmesh::ParallelMachine comm , std::istream & is )
{
  const unsigned p_size = phdmesh::parallel_machine_size( comm );
  std::string sval ;
  bool verify = false ;
  unsigned i = 2 ;
  unsigned j = 3 ;
  unsigned k = 1 ;

  if ( is.good() ) {
    is >> i ; is >> j ; is >> k ;
    is >> sval ;
  }

  if ( is.good() && sval == std::string("verify") ) {
    verify = true ;
    is >> sval ;
  }

  std::ostringstream exo_file_name ;

  if ( is.good() && sval == std::string("output") ) {
    exo_file_name << sval << "_np" << p_size << ".exo" ;
  }

  test_gears( comm , i , j , k , exo_file_name.str() , verify );
}

