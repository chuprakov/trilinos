/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/TopologicalMetaData.hpp>

#include <stk_rebalance/Rebalance.hpp>
#include <stk_rebalance/Partition.hpp>
#include <stk_rebalance/ZoltanPartition.hpp>

typedef stk::mesh::Field<double> ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> VectorField ;

enum { nx = 2, ny = 2 };


STKUNIT_UNIT_TEST(UnitTestZoltanSimple, testUnit)
{
#ifdef STK_HAS_MPI
  stk::ParallelMachine comm(MPI_COMM_WORLD);
#else
  stk::ParallelMachine comm(0);
#endif

  unsigned spatial_dimension = 2;
  stk::mesh::MetaData meta_data( stk::mesh::TopologicalMetaData::entity_rank_names(spatial_dimension) );
  stk::mesh::BulkData bulk_data( meta_data , comm , 100 );
  stk::mesh::TopologicalMetaData top_data( meta_data, spatial_dimension );
  stk::mesh::Part & quad_part( top_data.declare_part<shards::Quadrilateral<4> >( "quad" ) );
  VectorField & coord_field( meta_data.declare_field< VectorField >( "coordinates" ) );
  ScalarField & weight_field( meta_data.declare_field< ScalarField >( "element_weights" ) );

  stk::mesh::put_field( coord_field , top_data.node_rank , meta_data.universal_part() );
  stk::mesh::put_field(weight_field , top_data.element_rank , meta_data.universal_part() );

  meta_data.commit();

  //const unsigned p_size = bulk_data.parallel_size();
  const unsigned p_rank = bulk_data.parallel_rank();

  bulk_data.modification_begin();

  if ( p_rank == 0 ) {
    const unsigned nnx = nx + 1 ;
    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::EntityId nodes[4] ;
        nodes[0] = 1 + ix + iy * nnx ;
        nodes[1] = 2 + ix + iy * nnx ;
        nodes[2] = 2 + ix + ( iy + 1 ) * nnx ;
        nodes[3] = 1 + ix + ( iy + 1 ) * nnx ;

        stk::mesh::declare_element( bulk_data , quad_part , elem , nodes );
      }
    }

    for ( unsigned iy = 0 ; iy < ny ; ++iy ) {
      for ( unsigned ix = 0 ; ix < nx ; ++ix ) {
        stk::mesh::EntityId elem = 1 + ix + iy * nx ;
        stk::mesh::Entity * e = bulk_data.get_entity( top_data.element_rank, elem );
        double * const e_weight = stk::mesh::field_data( weight_field , *e );
        *e_weight = 1.0;
      }
    }
  }

  // Only P0 has any nodes or elements
  if ( p_rank == 0 ) {
    STKUNIT_ASSERT( ! bulk_data.buckets( top_data.node_rank ).empty() );
    STKUNIT_ASSERT( ! bulk_data.buckets( top_data.element_rank ).empty() );
  }
  else {
    STKUNIT_ASSERT( bulk_data.buckets( top_data.node_rank ).empty() );
    STKUNIT_ASSERT( bulk_data.buckets( top_data.element_rank ).empty() );
  }

  bulk_data.modification_end();

  Teuchos::ParameterList emptyList;
  stk::rebalance::Zoltan zoltan_partition(comm, spatial_dimension, emptyList);

  stk::mesh::Selector selector(meta_data.universal_part());

  // Force a rebalance by using imblance_threshhold < 1.0
  double imblance_threshhold = 0.5;
  bool do_rebal = stk::rebalance::rebalance_needed(bulk_data, meta_data, weight_field, comm, imblance_threshhold);
  if( do_rebal )
    stk::rebalance::rebalance(bulk_data, selector, NULL, &weight_field, zoltan_partition);
}
