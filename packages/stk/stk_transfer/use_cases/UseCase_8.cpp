/*------------------------------------------------------------------------*/
/*                 Copyright 2013 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <Intrepid_FieldContainer.hpp>
#include <boost/shared_ptr.hpp>


#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinMesh.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <stk_io/StkMeshIoBroker.hpp>
#include <init/Ionit_Initializer.h>

#include <stk_util/use_cases/UseCaseEnvironment.hpp>
#include <stk_util/diag/PrintTimer.hpp>

#include "STKElem.hpp"
#include "MDMesh.hpp"
#include "FEInterpolation.hpp"
#include <stk_transfer/GeometricTransfer.hpp>

namespace bopt = boost::program_options;

typedef stk::mesh::Field<double>                       ScalarField ;
typedef stk::mesh::Field<double, stk::mesh::Cartesian> CartesianField ;


bool use_case_8_driver(stk::ParallelMachine  comm,
                      const std::string &working_directory,
                      const std::string &domain_mesh,
                      const std::string &domain_filetype)
{
  stk::diag::Timer timer("Transfer Use Case 8",
                          use_case::TIMER_TRANSFER,
                          use_case::timer());
  stk::diag::Timer timer_node_to_node(" Element To Point", timer);
  use_case::timerSet().setEnabledTimerMask(use_case::TIMER_ALL);

  bool status = true;

  enum {           DIM = 3  };
  const double TOLERANCE = 0.000001;
  const double  rand_max = RAND_MAX;
  enum {   TONUMPOINTS = 100  };

  typedef Intrepid::FieldContainer<double>  MDArray;

  MDArray ToPoints   (  TONUMPOINTS,DIM),
          ToValues   (  TONUMPOINTS,  1);
  for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
    for (unsigned j=0 ; j<DIM; ++j) {
      ToPoints(i,j) = rand()/rand_max;
    }
  }

  const stk::mesh::EntityRank node_rank = stk::mesh::MetaData::NODE_RANK;
  const stk::mesh::EntityRank elem_rank = stk::mesh::MetaData::ELEMENT_RANK;

  const std::string data_field_name = "Sum_Of_Coordinates";

  stk::io::StkMeshIoBroker domain_mesh_data(comm);
  const std::string filename = working_directory + domain_mesh;
  domain_mesh_data.open_mesh_database(filename, domain_filetype);
  domain_mesh_data.create_input_mesh();

  stk::mesh::MetaData &domain_meta_data = domain_mesh_data.meta_data();
  stk::mesh::Part & domain_block        = domain_meta_data.declare_part("elements", elem_rank);
  stk::mesh::CellTopology hex_top (shards::getCellTopologyData<shards::Hexahedron<> >());
  stk::mesh::CellTopology quad_top(shards::getCellTopologyData<shards::Quadrilateral<> >());
  const stk::mesh::EntityRank side_rank    = domain_meta_data.side_rank();
  stk::mesh::Part & block_skin       = domain_meta_data.declare_part("skin", side_rank);

  stk::mesh::set_cell_topology( domain_block,  hex_top );
  stk::mesh::set_cell_topology( block_skin,    quad_top );

  ScalarField &domain_coord_sum_field = stk::mesh::put_field(
                        domain_meta_data.declare_field<ScalarField>(data_field_name),
                        node_rank ,
                        domain_meta_data.universal_part() );
  domain_meta_data.commit();

  domain_mesh_data.populate_bulk_data();
  stk::mesh::BulkData &domain_bulk_data = domain_mesh_data.bulk_data();
  stk::mesh::PartVector add_parts(1,&block_skin);
  stk::mesh::skin_mesh(domain_bulk_data, add_parts);
  // For this use case, the domain consists of an axis-aligned
  // bounding box for each 'domain_entity' in the mesh.  The range is a
  // PointBoundingBox3D at the centroid of each 'range_entity'.  The id of the point
  // will be the same as the id of the containing entity.  If the
  // mesh contains solid elements only, and the range_mesh matches the
  // domain_mesh, then the search should return a single box for each
  // point and the id of the box should match the id of the point.

  CartesianField const& domain_coord_field = static_cast<CartesianField const&>(domain_mesh_data.get_coordinate_field());

  stk::mesh::Selector locally_owned= domain_meta_data.locally_owned_part();

  {
    std::vector<stk::mesh::Entity> node_entities;
    stk::mesh::get_selected_entities(locally_owned, domain_bulk_data.buckets(node_rank), node_entities);
    const size_t num_entities = node_entities.size();
    for (size_t i = 0; i < num_entities; ++i) {
      const stk::mesh::Entity entity = node_entities[i];
      double *entity_coordinates = domain_bulk_data.field_data(domain_coord_field, entity);
      double *entity_coord_sum   = domain_bulk_data.field_data(domain_coord_sum_field, entity);
      *entity_coord_sum = entity_coordinates[0] + entity_coordinates[1] + entity_coordinates[2];
    }
  }

  std::vector<stk::mesh::Entity> domain_entities;
  stk::mesh::get_selected_entities(locally_owned, domain_bulk_data.buckets(elem_rank), domain_entities);

  const double radius=.25;
  const std::vector<stk::mesh::FieldBase*> from_fields(1, &domain_coord_sum_field);
  boost::shared_ptr<stk::transfer::STKElem >
    transfer_domain_mesh (new stk::transfer::STKElem(domain_entities, domain_coord_field, from_fields));
  boost::shared_ptr<stk::transfer:: MDMesh >
    transfer_range_mesh  (new stk::transfer:: MDMesh(ToValues, ToPoints,   radius, comm));


  stk::transfer::GeometricTransfer<
    class stk::transfer::FEInterpolate<
      class stk::transfer::STKElem,
      class stk::transfer::MDMesh
    >
  >
  transfer(transfer_domain_mesh, transfer_range_mesh, "STK Transfer test Use case 8");

  {
    stk::diag::TimeBlock __timer_node_to_node(timer_node_to_node);
    try {
      transfer.initialize();
      transfer.apply();
    } catch (std::exception &e) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an std::exception with what string:"
                <<e.what()
                <<"      rethrowing....."
                <<std::endl;
      status = status && false;
    } catch (...) {
      std::cout <<__FILE__<<":"<<__LINE__
                <<" Caught an exception, rethrowing..."
                <<std::endl;
      status = status && false;
    }
  }

  if (status) {

    bool success = true;
    for (unsigned i=0 ; i<TONUMPOINTS; ++i) {
      double check_l = 0;
      for (unsigned j=0 ; j<DIM; ++j) check_l += ToPoints(i,j);
      if (TOLERANCE < fabs(check_l-ToValues(i,0))) {
        std::cout <<__FILE__<<":"<<__LINE__
                  <<" EntityKey:"<<i
                  <<" ToPoints:"<<ToPoints(i,0)<<" "<<ToPoints(i,1)<<" "<<ToPoints(i,2)
                  <<" ToValues:"<<ToValues(i,0)
                  <<" check:"<<check_l
                  <<" error:"<<fabs(check_l-ToValues(i,0))
                  <<std::endl;
        success = false;
      }
    }
    status = status && success;
  }
  timer.stop();
//stk::diag::printTimersTable(std::cout, timer,
//      stk::diag::METRICS_CPU_TIME | stk::diag::METRICS_WALL_TIME, false, comm);


  const bool collective_result = use_case::print_status(comm, status);
  return collective_result;
}
