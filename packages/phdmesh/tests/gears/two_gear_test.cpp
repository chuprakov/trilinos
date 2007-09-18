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

#include "Gears.hpp"

using namespace phdmesh ;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void test_two_gears( ParallelMachine pm , std::istream & )
{
  const unsigned p_rank = parallel_machine_rank( pm );

  const unsigned kernel_capacity[ EntityTypeMaximum ] =
    { 100 , 100 , 100 , 100 , 100 };
    // { 20 , 20 , 20 , 20 , 20 };

  if ( p_rank == 0 ) {
    std::cout << "sizeof(unsigned long) = "
              <<  sizeof(unsigned long) << std::endl ;
  }

  Schema S( 3 , pm);

  GearFields gear_fields( S );

  //------------------------------
  // Two circular gears with a 'z' axis of rotation.
  // They touch as the origin through the thickness

  const double center_A[3] = {  1.0 , 0.0 , 0.0 };
  const double center_B[3] = { -1.0 , 0.0 , 0.0 };
  const double rad_min = 0.7 ;
  const double rad_max = 1.0 ;
  const double z_min   = -0.1 ;
  const double z_max   =  0.1 ;
  const unsigned z_num   = 3 ; // Number node planes through thickness
  const unsigned rad_num = 4 ; // Number of node shells through radius

#if 0
  enum { NAngle_P = 12 };
  const unsigned p_size    = parallel_machine_size( pm );
  const unsigned angle_num = NAngle_P * p_size ; // Number around radius
#else
  const unsigned angle_num = 72 ;
#endif
  
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

  //------------------------------

  S.commit();

  Mesh M(S,kernel_capacity);

  A.mesh( M );
  B.mesh( M );

  comm_mesh_discover_sharing( M );
  comm_mesh_regenerate_aura( M );

  if ( ! comm_mesh_verify_parallel_consistency( M ) ) {
    return ;
  }

  {
    unsigned long counts[ EntityTypeMaximum ];
    unsigned long max_id[ EntityTypeMaximum ];

    comm_mesh_stats( M , counts , max_id );

    if ( p_rank == 0 ) {
      std::cout << "TWO_GEARS Global Counts = " 
                << "node = " << counts[0] << " "
                << "edge = " << counts[1] << " "
                << "face = " << counts[2] << " "
                << "elem = " << counts[3] << " "
                << "other = " << counts[4] << std::endl ;
      std::cout << "TWO_GEARS Global MaxId  = " 
                << "node = " << max_id[0] << " "
                << "edge = " << max_id[1] << " "
                << "face = " << max_id[2] << " "
                << "elem = " << max_id[3] << " "
                << "other = " << max_id[4] << std::endl ;
    }

    partset_entity_count( M , S.owns_part() , counts );

    std::cout << "TWO_GEARS P" << p_rank << " Owned Counts = " 
              << "node = " << counts[0] << " "
              << "edge = " << counts[1] << " "
              << "face = " << counts[2] << " "
              << "elem = " << counts[3] << " "
              << "other = " << counts[4] << std::endl ;

    partset_entity_count( M , S.shares_part() , counts );

    std::cout << "TWO_GEARS P" << p_rank << " Shared Counts = " 
              << "node = " << counts[0] << " "
              << "edge = " << counts[1] << " "
              << "face = " << counts[2] << " "
              << "elem = " << counts[3] << " "
              << "other = " << counts[4] << std::endl ;

    partset_entity_count( M , S.aura_part() , counts );

    std::cout << "TWO_GEARS P" << p_rank << " Aura Counts = " 
              << "node = " << counts[0] << " "
              << "edge = " << counts[1] << " "
              << "face = " << counts[2] << " "
              << "elem = " << counts[3] << " "
              << "other = " << counts[4] << std::endl ;
  }

  //------------------------------

  if ( p_rank == 0 ) {
    std::cout << "TWO_GEARS test successful" << std::endl ;
    std::cout.flush();
  }
}


