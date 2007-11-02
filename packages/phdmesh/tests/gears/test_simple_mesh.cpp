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

#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>

#include <util/ParallelComm.hpp>

#include <mesh/Schema.hpp>
#include <mesh/Mesh.hpp>
#include <mesh/Comm.hpp>

using namespace phdmesh ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_simple_mesh( ParallelMachine pm , std::istream & )
{
  static const char method[] = "test_simple_mesh" ;

  const unsigned p_rank = parallel_machine_rank( pm );
  const unsigned p_size = parallel_machine_size( pm );

  //--------------------------------------------------------------------
  // Define a mesh schema: the parts and fields.

  Schema S ;

  // Nodal coordinate field dimensioned to 3 everywhere in the mesh

  Field<double,1> & node_coordinates =
    S.declare_field<double,1>( Node , std::string("coordinates") );

  S.declare_field_dimension( node_coordinates , S.universal_part() , 3 );
  
  // Done defining the schema, commit it.

  S.commit();

  //--------------------------------------------------------------------
  // Create mesh bulk data conformal to the schema.

  // Maximum number of entries in a mesh kernel,
  // { nodes , edges, faces , elements , other }

  const unsigned kernel_capacity[ EntityTypeMaximum ] =
    { 100 , 100 , 100 , 100 , 100 };

  Mesh M( S , pm , kernel_capacity );

  // Define a trivial mesh, one hex element per processor
  // Node and element identifiers must not be zero.

  // Base of hex
  const unsigned long node_id_1 = p_rank * 4 + 1 ;
  const unsigned long node_id_2 = node_id_1 + 1 ;
  const unsigned long node_id_3 = node_id_1 + 2 ;
  const unsigned long node_id_4 = node_id_1 + 3 ;

  // Top of hex
  const unsigned long node_id_5 = node_id_1 + 4 ;
  const unsigned long node_id_6 = node_id_1 + 5 ;
  const unsigned long node_id_7 = node_id_1 + 6 ;
  const unsigned long node_id_8 = node_id_1 + 7 ;

  const unsigned long elem_id = p_rank + 1 ;

  std::vector<Part*> parts ;
  { // Assume this processor owns everything it declares,
    // will resolve parallel sharing later.
    Part * const owns_part = & S.owns_part();
    parts.push_back(owns_part);
  }

  // Declare node and element entities:

  Entity & node_1 = M.declare_entity( Node , node_id_1 , parts , p_rank );
  Entity & node_2 = M.declare_entity( Node , node_id_2 , parts , p_rank );
  Entity & node_3 = M.declare_entity( Node , node_id_3 , parts , p_rank );
  Entity & node_4 = M.declare_entity( Node , node_id_4 , parts , p_rank );
  Entity & node_5 = M.declare_entity( Node , node_id_5 , parts , p_rank );
  Entity & node_6 = M.declare_entity( Node , node_id_6 , parts , p_rank );
  Entity & node_7 = M.declare_entity( Node , node_id_7 , parts , p_rank );
  Entity & node_8 = M.declare_entity( Node , node_id_8 , parts , p_rank );

  Entity & elem = M.declare_entity( Element , elem_id , parts , p_rank );

  // Declare element <-> node connections,
  // these are required to be unique by 'method'

  M.declare_connection( elem , node_1 , 1 , method );
  M.declare_connection( elem , node_2 , 2 , method );
  M.declare_connection( elem , node_3 , 3 , method );
  M.declare_connection( elem , node_4 , 4 , method );
  M.declare_connection( elem , node_5 , 5 , method );
  M.declare_connection( elem , node_6 , 6 , method );
  M.declare_connection( elem , node_7 , 7 , method );
  M.declare_connection( elem , node_8 , 8 , method );

  // Set node coordinates:

  double * const node_1_coord = node_1.data( node_coordinates );
  double * const node_2_coord = node_2.data( node_coordinates );
  double * const node_3_coord = node_3.data( node_coordinates );
  double * const node_4_coord = node_4.data( node_coordinates );
  double * const node_5_coord = node_5.data( node_coordinates );
  double * const node_6_coord = node_6.data( node_coordinates );
  double * const node_7_coord = node_7.data( node_coordinates );
  double * const node_8_coord = node_8.data( node_coordinates );

  node_1_coord[0] = 0 ; node_1_coord[1] = 0 ; node_1_coord[2] = p_rank ;
  node_2_coord[0] = 1 ; node_2_coord[1] = 0 ; node_2_coord[2] = p_rank ;
  node_3_coord[0] = 1 ; node_3_coord[1] = 1 ; node_3_coord[2] = p_rank ;
  node_4_coord[0] = 0 ; node_4_coord[1] = 1 ; node_4_coord[2] = p_rank ;
  node_5_coord[0] = 0 ; node_5_coord[1] = 0 ; node_5_coord[2] = p_rank + 1 ;
  node_6_coord[0] = 1 ; node_6_coord[1] = 0 ; node_6_coord[2] = p_rank + 1 ;
  node_7_coord[0] = 1 ; node_7_coord[1] = 1 ; node_7_coord[2] = p_rank + 1 ;
  node_8_coord[0] = 0 ; node_8_coord[1] = 1 ; node_8_coord[2] = p_rank + 1 ;

  // Determine proper parallel sharing and ownership
  comm_mesh_discover_sharing( M );

  // Generate the parallel ghosting 'aura'
  comm_mesh_regenerate_aura( M );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << method << " is not parallel consistent" << std::endl ;
    }
    return ;
  }

  // Get the global counts and identifier stats
  {
    unsigned long counts[ EntityTypeMaximum ];
    unsigned long max_id[ EntityTypeMaximum ];

    comm_mesh_stats( M , counts , max_id );

    if ( p_rank == 0 ) {
      std::cout << method
                << " Stats for { node , edge , face , element , other }"
                << std::endl ;
      std::cout << "  Global Counts = {" 
                << " " << counts[0]
                << " " << counts[1]
                << " " << counts[2]
                << " " << counts[3]
                << " " << counts[4]
                << " }" << std::endl ;
      std::cout << "  Global MaxId  = {" 
                << " " << max_id[0]
                << " " << max_id[1]
                << " " << max_id[2]
                << " " << max_id[3]
                << " " << max_id[4]
                << " }" << std::endl ;
    }

    for ( unsigned p = 0 ; p < p_size ; ++p ) {
      parallel_machine_barrier( pm );
      if ( p_rank == p ) {
        partset_entity_count( M , S.owns_part() , counts );

        std::cout << "  P" << p_rank << " Owned  Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " }" << std::endl ;

        partset_entity_count( M , S.shares_part() , counts );

        std::cout << "  P" << p_rank << " Shared Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " }" << std::endl ;

        partset_entity_count( M , S.aura_part() , counts );

        std::cout << "  P" << p_rank << " Aura   Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " }" << std::endl ;

        std::cout.flush();
      }
      parallel_machine_barrier( pm );
    }
  }

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << method << " successful" << std::endl ;
    std::cout.flush();
  }
}


