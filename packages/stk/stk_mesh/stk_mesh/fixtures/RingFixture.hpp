/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_FIXTURES_RING_FIXTURE_HPP
#define STK_MESH_FIXTURES_RING_FIXTURE_HPP

#include <stk_util/parallel/Parallel.hpp>

#include <stk_util/environment/ReportHandler.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <Shards_BasicTopologies.hpp>

namespace stk {
namespace mesh {
namespace fixtures {

/**
 * Creates a ring mesh (circular loop of elements and nodes). Note that we create
 * a part for each locally owned element. This fixture is 1d, so elements are just lines.
 */

class RingFixture {
 public:
  const int             m_spatial_dimension;
  MetaData              m_meta_data;
  BulkData              m_bulk_data;
  PartVector            m_element_parts ;
  Part &                m_element_part_extra ;
  const size_t          m_num_element_per_proc ;
  std::vector<EntityId> m_node_ids , m_element_ids ;

  RingFixture( stk::ParallelMachine pm ,
               unsigned num_element_per_proc = 10 ,
               bool use_element_parts = false );

  ~RingFixture() {}

  /**
   * Generate a simple loop of mesh entities:
   * node[i] : element[i] : node[ ( i + 1 ) % node.size() ]
   */
  void generate_mesh();

  /**
   * Make sure that element->owner_rank() == element->node[1]->owner_rank()
   */
  void fixup_node_ownership();

 private:

   RingFixture();
   RingFixture( const RingFixture & );
   RingFixture & operator = ( const RingFixture & );
};

}
}
}

#endif
