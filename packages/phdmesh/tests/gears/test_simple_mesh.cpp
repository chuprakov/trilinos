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

  // Get some of the predefined parts for later use...
  Part * const owns_part = & S.owns_part();
  Part * const univ_part = & S.universal_part();

  // Declare a part for the element block and side set,
  // these are automatically a subset of the universal part.

  Part * const elem_part = & S.declare_part( std::string("element_block") );
  Part * const face_part = & S.declare_part( std::string("side_set") );

  // Nodal coordinate field dimensioned to 3 everywhere in the mesh

  Field<double,1> & node_coordinates =
    S.declare_field<double,1>( Node , std::string("coordinates") );

  S.declare_field_dimension( node_coordinates , *univ_part , 3 );
  
  // Done defining the schema, commit it.

  S.commit();

  //--------------------------------------------------------------------
  // Create mesh bulk data conformal to the schema.

  // Maximum number of entries in a mesh kernel,
  // { nodes , edges, faces , elements , other }

  const unsigned kernel_capacity[ EntityTypeMaximum ] =
    { 100 , 100 , 100 , 100 , 100 };

  Mesh M( S , pm , kernel_capacity );

  // Define a trivial mesh, stack of hex elements
  // with one hex element per processor ordered by
  // processor rank.
  // Attach the +X face, use ExodusII element-node ordering
  // and element-face orientation.
  // Node and element identifiers must not be zero,
  // an identifier of zero is reserved for 'undefined'.

  // Base of this processor's hex
  const unsigned long node_id_1 = p_rank * 4 + 1 ;
  const unsigned long node_id_2 = node_id_1 + 1 ;
  const unsigned long node_id_3 = node_id_1 + 2 ;
  const unsigned long node_id_4 = node_id_1 + 3 ;

  // Top of this processor's hex
  const unsigned long node_id_5 = node_id_1 + 4 ;
  const unsigned long node_id_6 = node_id_1 + 5 ;
  const unsigned long node_id_7 = node_id_1 + 6 ;
  const unsigned long node_id_8 = node_id_1 + 7 ;

  const unsigned long elem_id = p_rank + 1 ;
  const unsigned long face_id = p_rank + 1 ;

  // Part membership for the elements and nodes:
  // 'owns_part'  Assume this processor owns everything it declares,
  //              will resolve parallel sharing later.

  std::vector<Part*> add_parts , remove_parts ;
  add_parts.push_back( owns_part );
  add_parts.push_back( elem_part );

  // Declare node and element entities:

  Entity & elem = M.declare_entity( Element , elem_id , add_parts , p_rank );

  Entity & node_1 = M.declare_entity( Node , node_id_1 , add_parts , p_rank );
  Entity & node_2 = M.declare_entity( Node , node_id_2 , add_parts , p_rank );
  Entity & node_3 = M.declare_entity( Node , node_id_3 , add_parts , p_rank );
  Entity & node_4 = M.declare_entity( Node , node_id_4 , add_parts , p_rank );
  Entity & node_5 = M.declare_entity( Node , node_id_5 , add_parts , p_rank );
  Entity & node_6 = M.declare_entity( Node , node_id_6 , add_parts , p_rank );
  Entity & node_7 = M.declare_entity( Node , node_id_7 , add_parts , p_rank );
  Entity & node_8 = M.declare_entity( Node , node_id_8 , add_parts , p_rank );

  // Declare element <-> node connections
  // These are required to have unique identifiers
  // by providing the 'method' argument.
  // If non-unique then an exception is thrown that includes
  // the text contained in the 'method' string.

  M.declare_connection( elem , node_1 , 1 , method );
  M.declare_connection( elem , node_2 , 2 , method );
  M.declare_connection( elem , node_3 , 3 , method );
  M.declare_connection( elem , node_4 , 4 , method );
  M.declare_connection( elem , node_5 , 5 , method );
  M.declare_connection( elem , node_6 , 6 , method );
  M.declare_connection( elem , node_7 , 7 , method );
  M.declare_connection( elem , node_8 , 8 , method );

  // Declare the face entity:

  add_parts.push_back( face_part );
  
  Entity & face = M.declare_entity( Face , face_id , add_parts , p_rank );

  // Declare element <-> face connection

  M.declare_connection( elem , face , 2 , method );

  // Declare face <-> node connections

  M.declare_connection( face , node_1 , 1 , method );
  M.declare_connection( face , node_2 , 2 , method );
  M.declare_connection( face , node_6 , 3 , method );
  M.declare_connection( face , node_5 , 4 , method );

  // Update the nodes on the face to also be members of the face part.

  M.change_entity_parts( node_1 , add_parts , remove_parts );
  M.change_entity_parts( node_2 , add_parts , remove_parts );
  M.change_entity_parts( node_6 , add_parts , remove_parts );
  M.change_entity_parts( node_5 , add_parts , remove_parts );

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

  // Verify that the parallel sharing and aura were generated properly
  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    if ( p_rank == 0 ) {
      std::cout << method
                << " FAILED: is not parallel consistent"
                << std::endl ;
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
        partset_entity_count( M , S.uses_part() , counts );

        std::cout << "  P" << p_rank << " Uses  Counts = {" 
                  << " " << counts[0]
                  << " " << counts[1]
                  << " " << counts[2]
                  << " " << counts[3]
                  << " " << counts[4]
                  << " }" << std::endl ;

        partset_entity_count( M , S.owns_part() , counts );

        std::cout << "  P" << p_rank << " Owns  Counts = {" 
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
  parallel_machine_barrier( pm );

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << method << " successful" << std::endl ;
    std::cout.flush();
  }
}


