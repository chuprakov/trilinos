/*------------------------------------------------------------------------*/
/*                 Copyright 2014 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stk_mesh/fixtures/degenerate_mesh.hpp>
#include <Shards_BasicTopologies.hpp>   // for Hexahedron
#include <sstream>                      // for ostringstream, etc
#include <stdexcept>                    // for runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/FEMHelpers.hpp>  // for declare_element
#include <stk_mesh/base/FindRestriction.hpp>  // for find_restriction
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include "Shards_CellTopologyTraits.hpp"
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase::Restriction, etc
#include "stk_mesh/base/Types.hpp"      // for EntityId
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }




namespace stk {
  namespace mesh {
    namespace fixtures {

      enum { SpatialDim = 3 };

      //----------------------------------------------------------------------

      //--------------------------------------------------------------------
      //

      void degenerate_mesh_meta_data(
				     stk::mesh::MetaData & meta_data ,
				     VectorFieldType & node_coord )
      {
	stk::mesh::Part & universal        = meta_data.universal_part();
	meta_data.declare_part_with_topology( "hexes", stk::topology::HEX_8);

	const stk::mesh::FieldBase::Restriction & res =
	  stk::mesh::find_restriction(node_coord, stk::topology::NODE_RANK , universal );

	if ( res.num_scalars_per_entity() != 3 ) {
	  std::ostringstream msg ;
	  msg << "stk_mesh/unit_tests/degenerate_mesh_meta_data FAILED, coordinate dimension must be 3 != "
	      << res.num_scalars_per_entity() ;
	  throw std::runtime_error( msg.str() );
	}
      }

      //--------------------------------------------------------------------
      /*----------------------------------------------------------------------
       * Degenerate hex mesh generation
       *
       * 12 hexes degenerated down to wedges.
       *
       *----------------------------------------------------------------------*/

      namespace {

	enum { node_count = 10 };
	enum { number_hex = 2 };

	/*  Z = 0 plane:
	 *
	 *    Y
	 *    ^ 
	 *    !
	 *     
	 *   4*       *5
	 *    |\     /|
	 *    | \   / |
	 *    |  \ /  |
	 *    *---*---* ----> X
	 *    1   2   3
	 *
	 *  Z = -1 plane:
	 *
	 *    Y
	 *    ^ 
	 *    !
	 *     
	 *   9*       *10
	 *    |\     /|
	 *    | \   / |
	 *    |  \ /  |
	 *    *---*---* ----> X
	 *    6   7   8
	 *
	 */

	static const double node_coord_data[ node_count ][ SpatialDim ] = {
	  { 0 , 0 , 0 } , { 1 , 0 , 0 } , { 2, 0, 0} , {0, 1, 0}, {2, 1, 0}, 
	  { 0 , 0 ,-1 } , { 1 , 0 ,-1 } , { 2, 0,-1} , {0, 1,-1}, {2, 1,-1}};

	static const stk::mesh::EntityId hex_node_ids[number_hex][ shards::Hexahedron<8> ::node_count ] = {
	  { 1, 2, 7, 6, 4, 2,  7, 9}, 
	  { 2, 3, 8, 7, 2, 5, 10, 7}};

      }
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------

      void degenerate_mesh_bulk_data(
				     stk::mesh::BulkData & bulk_data ,
				     const VectorFieldType & node_coord )
      {
	static const char method[] =
	  "stk_mesh::fixtures::heterogenous_mesh_bulk_data" ;

	bulk_data.modification_begin();

	const stk::mesh::MetaData & meta_data = stk::mesh::MetaData::get(bulk_data);

	stk::mesh::Part & hex_block        = * meta_data.get_part("hexes",method);

	unsigned elem_id = 1 ;

	for ( unsigned i = 0 ; i < number_hex ; ++i , ++elem_id ) {
	  stk::mesh::declare_element( bulk_data, hex_block, elem_id, hex_node_ids[i] );
	}

	for ( unsigned i = 0 ; i < node_count ; ++i ) {

	  stk::mesh::Entity const node = bulk_data.get_entity( stk::topology::NODE_RANK , i + 1 );

	  double * const coord = stk::mesh::field_data( node_coord , node );

	  coord[0] = node_coord_data[i][0] ;
	  coord[1] = node_coord_data[i][1] ;
	  coord[2] = node_coord_data[i][2] ;
	}

	bulk_data.modification_end();
      }

    }}}

//----------------------------------------------------------------------

