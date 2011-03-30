/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string.h>
#include <iostream>

#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>
#include <Ioss_NullEntity.h>

#include <stk_util/util/tokenize.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/FEMInterface.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <Shards_BasicTopologies.hpp>

namespace {

  const CellTopologyData *map_ioss_to_topology( const std::string &element_type ,
							const int node_count)
  {
    /// \todo REFACTOR Is there a good type to return for an "unknown"
    /// topology type other than NULL?

    const CellTopologyData* celltopo = NULL ;

    const char *etype = element_type.c_str();
    if ( 0 == strncasecmp( "circle" , etype , 6 ) ) {
      celltopo = shards::getCellTopologyData< shards::Particle >();
    }
    else if ( 0 == strncasecmp( "sphere" , etype , 6) ) {
      celltopo = shards::getCellTopologyData< shards::Particle >();
    }
    // bar, beam, truss, rod...
    else if ( 0 == strncasecmp( "bar" , etype , 3 ) ) {
      if ( node_count == 2 ) {
	celltopo = shards::getCellTopologyData< shards::Beam<2> >();
      }
      else if ( node_count == 3 ) {
	celltopo = shards::getCellTopologyData< shards::Beam<3> >();
      }
    }
    else if ( 0 == strncasecmp( "shellline2d" , etype , 11 ) ) {
      if ( node_count == 2) {
	celltopo = shards::getCellTopologyData< shards::ShellLine<2> >();
      }
      else if ( node_count == 3) {
	celltopo = shards::getCellTopologyData< shards::ShellLine<3> >();
      }
    } else if ( 0 == strncasecmp( "shell" , etype , 5 ) ) {
      // shell4, shell8, shell9
      if ( node_count == 4 ) {
	celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<4> >();
      }
      else if ( node_count == 8 ) {
	celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<8> >();
      }
      else if ( node_count == 9 ) {
	celltopo = shards::getCellTopologyData< shards::ShellQuadrilateral<9> >();
      }
    }
    else if ( 0 == strncasecmp( "quad" , etype , 3 ) ) {
      // The 2D types would be quad4, quad8, and quad9.
      // The 3D types would be quad faces of a hex... quadface4,
      // quadface8, quadface9.
      if ( node_count == 4 ) {
	celltopo = shards::getCellTopologyData< shards::Quadrilateral<4> >();
      }
      else if ( node_count == 8 ) {
	celltopo = shards::getCellTopologyData< shards::Quadrilateral<8> >();
      }
      else if ( node_count == 9 ) {
	celltopo = shards::getCellTopologyData< shards::Quadrilateral<9> >();
      }
    }
    else if ( 0 == strncasecmp( "trishell" , etype , 8 ) ) {
      if ( node_count == 3 ) {
	celltopo = shards::getCellTopologyData< shards::ShellTriangle<3> >();
      }
      else if ( node_count == 6 ) {
	celltopo = shards::getCellTopologyData< shards::ShellTriangle<6> >();
      }
    }

    else if (0 == strncasecmp("triface", etype, 7) ||
	     0 == strncasecmp("tri",     etype, 3)) {
      if ( node_count == 3 ) {
	celltopo = shards::getCellTopologyData< shards::Triangle<3> >();
      }
      else if ( node_count == 4 ) {
	celltopo = shards::getCellTopologyData< shards::Triangle<4> >();
      }
      else if ( node_count == 6 ) {
	celltopo = shards::getCellTopologyData< shards::Triangle<6> >();
      }
    }

    else if ( 0 == strncasecmp( "pyramid" , etype , 7 ) ) {
      if ( node_count == 5 ) {
	celltopo = shards::getCellTopologyData< shards::Pyramid<5> >();
      }
      else if ( node_count == 13 ) {
	celltopo = shards::getCellTopologyData< shards::Pyramid<13> >();
      }
      else if ( node_count == 14 ) {
	celltopo = shards::getCellTopologyData< shards::Pyramid<14> >();
      }
    }

    /// \todo REFACTOR Need to handle 8-node tet...
    else if ( 0 == strncasecmp( "tetra" , etype , 5 ) ) {
      if ( node_count == 4 ) {
	celltopo = shards::getCellTopologyData< shards::Tetrahedron<4> >();
      }
      else if ( node_count ==  8 ) {
	celltopo = shards::getCellTopologyData< shards::Tetrahedron<8> >();
      }
      else if ( node_count == 10 ) {
	celltopo = shards::getCellTopologyData< shards::Tetrahedron<10> >();
      }
    }
    else if ( 0 == strncasecmp( "wedge" , etype , 5 ) ) {
      if ( node_count == 6 ) {
	celltopo = shards::getCellTopologyData< shards::Wedge<6> >();
      }
      else if ( node_count == 15 ) {
	celltopo = shards::getCellTopologyData< shards::Wedge<15> >();
      }
      else if ( node_count == 18 ) {
	celltopo = shards::getCellTopologyData< shards::Wedge<18> >();
      }
    }
    else if ( 0 == strncasecmp( "hex" , etype , 3 ) ) {
      if ( node_count == 8 ) {
	celltopo = shards::getCellTopologyData< shards::Hexahedron<8> >();
      }
      else if ( node_count == 20 ) {
	celltopo = shards::getCellTopologyData< shards::Hexahedron<20> >();
      }
      else if ( node_count == 27 ) {
	celltopo = shards::getCellTopologyData< shards::Hexahedron<27> >();
      }
    }

    else if (0 == strncasecmp("edge", etype, 4)) {
      if ( node_count == 2) {
	// edge2, edge2d2, edge3d2
	celltopo = shards::getCellTopologyData< shards::Line<2> >();
      }
      else if ( node_count == 3) {
	// edge3, edge2d3, edge3d3
	celltopo = shards::getCellTopologyData< shards::Line<3> >();
      }
    }

    else if (0 == strncasecmp("node", etype, 4)) {
      celltopo = shards::getCellTopologyData< shards::Node >();
    }

    if ( NULL == celltopo ) {
      std::ostringstream oss;
      oss << "ERROR, unsupported topology name = '" << element_type
	  << "' , node_count = " << node_count;
      throw std::runtime_error(oss.str());
    }

    return celltopo;
  }

  const stk::mesh::FieldBase *declare_ioss_field(stk::mesh::MetaData &meta,
						 stk::mesh::EntityRank type,
						 stk::mesh::Part &part,
						 const Ioss::Field &io_field,
						 bool use_cartesian_for_scalar)
    {
      stk::mesh::FieldBase *field_ptr = NULL;
      std::string field_type = io_field.transformed_storage()->name();
      std::string name = io_field.get_name();
      int num_components = io_field.transformed_storage()->component_count();

      if (field_type == "scalar" || num_components == 1) {
	if (!use_cartesian_for_scalar) {
	  stk::mesh::Field<double> & field = meta.declare_field<stk::mesh::Field<double> >(name);
	  stk::mesh::put_field(field, type, part);
	  field_ptr = &field;
	} else {
	  stk::mesh::Field<double, stk::mesh::Cartesian> & field =
	    meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(name);
	  stk::mesh::put_field(field, type, part, 1);
	  field_ptr = &field;
	}
      }
      else if (field_type == "vector_2d") {
	stk::mesh::Field<double, stk::mesh::Cartesian> & field =
	  meta.declare_field<stk::mesh::Field<double, stk::mesh::Cartesian> >(name);
	stk::mesh::put_field(field, type, part, 2);
	field_ptr = &field;
      }
      else if (field_type == "vector_3d") {
	stk::mesh::Field<double, stk::mesh::Cartesian> & field =
	  meta.declare_field<stk::mesh::Field<double,
	  stk::mesh::Cartesian> >(name);
	stk::mesh::put_field(field, type, part, 3);
	field_ptr = &field;
      }
      else if (field_type == "sym_tensor_33") {
	stk::mesh::Field<double, stk::mesh::SymmetricTensor> & field =
	  meta.declare_field<stk::mesh::Field<double,
	  stk::mesh::SymmetricTensor> >(name);
	stk::mesh::put_field(field, type, part, 6);
	field_ptr = &field;
      }
      else if (field_type == "full_tensor_36") {
	stk::mesh::Field<double, stk::mesh::FullTensor> & field =
	  meta.declare_field<stk::mesh::Field<double,
	  stk::mesh::FullTensor> >(name);
	stk::mesh::put_field(field, type, part, 9);
	field_ptr = &field;
      }
      else {
	// Just create a field with the correct number of components...
	stk::mesh::Field<double,shards::ArrayDimension> & field =
	  meta.declare_field<stk::mesh::Field<double,shards::ArrayDimension> >(name);
	stk::mesh::put_field(field, type, part, num_components);
	field_ptr = &field;
      }

      if (field_ptr != NULL) {
	stk::io::set_field_role(*field_ptr, io_field.get_role());
      }
      return field_ptr;
    }

    template <typename T>
    void internal_field_data_from_ioss(const Ioss::Field &io_field,
				       const stk::mesh::FieldBase *field,
				       std::vector<stk::mesh::Entity*> &entities,
				       Ioss::GroupingEntity *io_entity,
				       T /*dummy */)
    {
      int field_component_count = io_field.transformed_storage()->component_count();

      std::vector<T> io_field_data;
      size_t io_entity_count = io_entity->get_field_data(io_field.get_name(), io_field_data);
      assert(io_field_data.size() == entities.size() * field_component_count);

      size_t entity_count = entities.size();

      if (io_entity_count != entity_count) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Field count mismatch for IO field '"
	       << io_field.get_name() << "'. The IO system has " << io_entity_count
	       << " entries, but the stk:mesh system has " << entity_count
	       << " entries. The two counts must match.";
	throw std::runtime_error(errmsg.str());
      }

      for (size_t i=0; i < entity_count; ++i) {
	/// \todo REFACTOR Is there a way to get the data from a
	/// "FieldBase*" field as a T* without the cast?
	if (entities[i] != NULL) {
	  T *fld_data = (T*)stk::mesh::field_data(*field, *entities[i]);
	  assert(fld_data != NULL);
	  for(int j=0; j<field_component_count; ++j) {
	    fld_data[j] = io_field_data[i*field_component_count+j];
	  }
	}
      }
    }

    template <typename T>
    void internal_field_data_to_ioss(const Ioss::Field &io_field,
				     const stk::mesh::FieldBase *field,
				     std::vector<stk::mesh::Entity*> &entities,
				     Ioss::GroupingEntity *io_entity,
				     T /*dummy */)
    {
      int field_component_count = io_field.transformed_storage()->component_count();
      size_t entity_count = entities.size();

      std::vector<T> io_field_data(entity_count*field_component_count);

      for (size_t i=0; i < entity_count; ++i) {
	/// \todo REFACTOR Is there a way to get the data from a
	/// "FieldBase*" field as a T* without the cast?
	if (entities[i] != NULL) {
	  T *fld_data = (T*)stk::mesh::field_data(*field, *entities[i]);
	  assert(fld_data != NULL);
	  for(int j=0; j<field_component_count; ++j) {
	    io_field_data[i*field_component_count+j] = fld_data[j];
	  }
	} else {
	  for(int j=0; j<field_component_count; ++j) {
	    io_field_data[i*field_component_count+j] = 0;
	  }
	}
      }

      size_t io_entity_count = io_entity->put_field_data(io_field.get_name(), io_field_data);
      assert(io_field_data.size() == entities.size() * field_component_count);

      if (io_entity_count != entity_count) {
	std::ostringstream errmsg;
	errmsg << "ERROR: Field count mismatch for IO field '"
	       << io_field.get_name() << "'. The IO system has " << io_entity_count
	       << " entries, but the stk:mesh system has " << entity_count
	       << " entries. The two counts must match.";
	throw std::runtime_error(errmsg.str());
      }
    }

}//namespace <empty>

namespace stk {
  namespace io {

    bool invalid_rank(stk::mesh::EntityRank rank)
    {
#ifdef USE_FEMMETADATA
      return rank == mesh::InvalidEntityRank;
#else
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
      return rank == (unsigned)mesh::EntityRankUndefined;
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
      return rank == mesh::InvalidEntityRank;
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
#endif
    }

    stk::mesh::EntityRank part_primary_entity_rank(const stk::mesh::Part &part)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(stk::mesh::MetaData::get(part).get_attribute<stk::mesh::fem::FEMMetaData>());

      if (mesh::MetaData::get(part).universal_part() == part) {
        if( fem_meta )
          {
            return stk::mesh::fem::NODE_RANK;
          }
        else
          {

#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
            return stk::mesh::Node;
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
            return stk::mesh::fem::NODE_RANK;
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */

          }
      }
      if (fem_meta)
        {
          return part.primary_entity_rank();
        }
      else
        {
#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return mesh::fem_entity_rank( part.primary_entity_rank() );
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          return part.primary_entity_rank();
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */

        }

    }

    stk::mesh::EntityRank element_rank(const stk::mesh::MetaData &meta)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        return fem_meta->element_rank();
      else
      {
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
        return stk::mesh::Element;
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta);
        return stk::mesh::fem::element_rank(fem);
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
      }

    }

    stk::mesh::EntityRank side_rank(const stk::mesh::MetaData &meta)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          // TODO return stk::mesh::fem::FEMMetaData::get(meta).side_rank();
          return stk::mesh::Face;
        }
      else
        {
#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return stk::mesh::Face;
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta);
          return stk::mesh::fem::side_rank(fem);
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        }
    }

    stk::mesh::EntityRank face_rank(const stk::mesh::MetaData &meta)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          // TODO
          //return stk::mesh::fem::FEMMetaData::get(meta).face_rank();

          return stk::mesh::Face;
        }
      else
        {
#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return stk::mesh::Face;
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta);
          return stk::mesh::fem::face_rank(fem);
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        }

    }

    stk::mesh::EntityRank edge_rank(const stk::mesh::MetaData &meta)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          // TODO
          // return stk::mesh::fem::FEMMetaData::get(meta).edge_rank();
          return stk::mesh::Edge;
        }
      else
        {
#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return stk::mesh::Edge;
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta);
          return stk::mesh::fem::edge_rank(fem);
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        }
    }

    stk::mesh::EntityRank node_rank(const stk::mesh::MetaData& meta)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          return stk::mesh::fem::FEMMetaData::get(meta).node_rank();
        }
      else
        {
#  ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return stk::mesh::Node;
#  else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          return stk::mesh::fem::NODE_RANK;
#  endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        }
    }

    void set_cell_topology(stk::mesh::Part &part, const CellTopologyData * const cell_topology)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(stk::mesh::MetaData::get(part).get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          stk::mesh::fem::set_cell_topology_new(part, cell_topology);
        }
      else
        {

#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          stk::mesh::set_cell_topology(part, cell_topology);
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          stk::mesh::fem::set_cell_topology(part, cell_topology);
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */

        }
    }

    const CellTopologyData *get_cell_topology(const stk::mesh::Part &part)
    {
      stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(stk::mesh::MetaData::get(part).get_attribute<stk::mesh::fem::FEMMetaData>());
      if( fem_meta )
        {
          return stk::mesh::fem::FEMMetaData::get(part).get_cell_topology(part).getCellTopologyData();
        }
      else
        {
#ifndef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
          return stk::mesh::get_cell_topology(part);
#else /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
          return stk::mesh::fem::get_cell_topology(part).getCellTopologyData();
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
        }
    }

    void initialize_spatial_dimension(stk::mesh::MetaData &meta, size_t spatial_dimension, const std::vector<std::string> &entity_rank_names)
    {
#ifdef SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS
      stk::mesh::fem::FEMInterface &fem = stk::mesh::fem::get_fem_interface(meta);

      meta.set_entity_rank_names(entity_rank_names);
      fem.set_spatial_dimension(spatial_dimension);
#endif /* SKIP_DEPRECATED_STK_MESH_TOPOLOGY_HELPERS */
    }

    void get_io_field_type(const stk::mesh::FieldBase *field, int num_comp, std::pair<std::string, Ioss::Field::BasicType> *result)
    {
      const std::string invalid("invalid");
      const std::string scalar("scalar");
      const std::string vector_2d("vector_2d");
      const std::string vector_3d("vector_3d");
      const std::string quaternion_2d("quaternion_2d");
      const std::string quaternion_3d("quaternion_3d");
      const std::string full_tensor_36("full_tensor_36");
      const std::string full_tensor_32("full_tensor_32");
      const std::string full_tensor_22("full_tensor_22");
      const std::string full_tensor_16("full_tensor_16");
      const std::string full_tensor_12("full_tensor_12");
      const std::string sym_tensor_33("sym_tensor_33");
      const std::string sym_tensor_31("sym_tensor_31");
      const std::string sym_tensor_21("sym_tensor_21");
      const std::string sym_tensor_13("sym_tensor_13");
      const std::string sym_tensor_11("sym_tensor_11");
      const std::string sym_tensor_10("sym_tensor_10");
      const std::string asym_tensor_03("asym_tensor_03");
      const std::string asym_tensor_02("asym_tensor_02");
      const std::string asym_tensor_01("asym_tensor_01");
      const std::string matrix_22("matrix_22");
      const std::string matrix_33("matrix_33");

      const unsigned rank = field->rank();
      const shards::ArrayDimTag * const * const tags = field->dimension_tags();

      result->second = Ioss::Field::INVALID;

      if ( field->type_is<double>() ) {
	result->second = Ioss::Field::REAL;
      }
      else if ( field->type_is<int>() ) {
	result->second = Ioss::Field::INTEGER;
      }

      if ( 0 == rank ) {
	result->first = scalar ;
      }
      else if ( 1 == rank ) {
	if ( tags[0] == & stk::mesh::Cartesian::tag() && 1 == num_comp ) {
	  result->first = scalar ;
	}
	else if ( tags[0] == & stk::mesh::Cartesian::tag() && 2 == num_comp ) {
	  result->first = vector_2d ;
	}
	else if ( tags[0] == & stk::mesh::Cartesian::tag() && 3 == num_comp ) {
	  result->first = vector_3d ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 9 == num_comp ) {
	  result->first = full_tensor_36 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 5 == num_comp ) {
	  result->first = full_tensor_32 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 4 == num_comp ) {
	  result->first = full_tensor_22 ;
	}
	else if ( tags[0] == & stk::mesh::FullTensor::tag() && 3 == num_comp ) {
	  result->first = full_tensor_12 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 6 == num_comp ) {
	  result->first = sym_tensor_33 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 4 == num_comp ) {
	  result->first = sym_tensor_31 ;
	}
	else if ( tags[0] == & stk::mesh::SymmetricTensor::tag() && 3 == num_comp ) {
	  result->first = sym_tensor_21 ;
	}
      }

      if ( result->first.empty() ) {
	std::ostringstream tmp ;
	tmp << "Real[" << num_comp << "]" ;
	result->first = tmp.str();
      }
    }

    //----------------------------------------------------------------------
    const Ioss::GroupingEntity *get_associated_ioss_entity(const mesh::Part &part)
    {
      const Ioss::GroupingEntity *entity = part.attribute<Ioss::GroupingEntity>();
      if (!entity || entity->type() == Ioss::INVALID_TYPE) {
	return NULL;
      } else {
	return entity;
      }
    }

    void put_io_part_attribute(mesh::Part & part, Ioss::GroupingEntity *entity)
    {
      if (part.attribute<Ioss::GroupingEntity>() != NULL) {
	std::string msg = "stk::io::put_io_part_attribute( ";
	msg += part.name();
	msg += " ) FAILED:";
	msg += " io_part_attribute is already defined";
	throw std::runtime_error( msg );
      }

      mesh::MetaData & meta = mesh::MetaData::get(part);
      if (entity) {
	meta.declare_attribute_no_delete(part, entity);
      } else {
	Ioss::GroupingEntity *attr = new Ioss::NullEntity();
	meta.declare_attribute_with_delete(part, attr);
      }
    }

    /** Determine whether the field is defined on the specified part
     * and should also be associated with a Ioss GroupingEntity for
     * input or output
     */
    bool is_valid_part_field(const stk::mesh::FieldBase *field,
			     stk::mesh::EntityRank part_type,
			     stk::mesh::Part &part,
			     stk::mesh::Part &universal,
			     Ioss::Field::RoleType filter_role,
			     bool add_all)
    {
      const Ioss::Field::RoleType *role = stk::io::get_field_role(*field);
      if (!add_all && (role == NULL || *role != filter_role)) {
	return false;
      }

      const stk::mesh::FieldBase::Restriction &res = field->restriction(part_type, part);
      if (res.stride[0] > 0) {
	// The field exists on the current 'part'.  Now check (for
	// node types only) whether the 'part' is *either* the
	// universal_part() *or* the field *doesn't* exist on the
	// universal part...
	// Note that for "node" type parts, the IO database has a part
	// (nodeblock) that corresponds to the universal_part(), so
	// fields defined on all nodes are output on the nodeblock and
	// fields defined on only a subset of the parts should be
	// output on that subset which maps to a nodeset.  For other
	// part types, there is no IO entity that corresponds to a
	// universal part, so fields must be output on the part they
	// exist on.  There may be a problem if we start using element
	// sets ..., but wait until we get to that point; current code
	// works with current entity set.
	if (part_type != node_rank(mesh::MetaData::get(part)) || part == universal) {
	  return true;
	}

	const stk::mesh::FieldBase::Restriction &res_universe = field->restriction(part_type, universal);
	if (res_universe.stride[0] <= 0) {
	  // Field exists on current part, but not on the universal
	  // set (and this part is not the universal part)
	  return true;
	}
      }
      return false;
    }

    /// \todo QUESTION Should this function be at the application level,
    /// or provided by stk_io? In either case, applications should have
    /// capabilty to register new mappings.
    // ========================================================================
    const CellTopologyData *map_topology_ioss_to_cell(const Ioss::ElementTopology *topology)
    {
      /// \todo REFACTOR Consider either using or integrating the
      /// Trilinos CellTopology package into or with the
      /// Ioss::ElementTopology classes. That would then totally
      /// eliminate the need for these fragile mapping functions.
      /// However, it would still need to be extensible via application
      /// registration of new type mappings.

      std::string name         = topology->name();
      int io_nodes_per_element = topology->number_nodes();

      const CellTopologyData *cell_topology = map_ioss_to_topology(name, io_nodes_per_element);

      return cell_topology;
    }

    std::string map_topology_cell_to_ioss( const CellTopologyData *cell_top,
					   int spatial_dimension)
    {
      std::string extype = "unknown";

      if (cell_top == NULL)
	return extype;

      if(strcmp(cell_top->name, "super") == 0) {
          std::stringstream oss;
          oss << "super" << cell_top->node_count;
          return oss.str();
      }

      switch( cell_top->key ) {
      case shards::Node::key : extype.assign( "node" ); break ;
      case shards::Particle::key :
	if (spatial_dimension == 2) extype = "circle1";
	else                        extype = "sphere1";
	break ;

      case shards::Beam<2>::key :
	extype = "bar2";
	break ;

      case shards::Beam<3>::key :
	extype = "bar3";
	break ;

      case shards::Line<2>::key :
	if (spatial_dimension == 2) extype = "edge2d2";
	else                        extype = "edge2";
	break ;
      case shards::Line<3>::key :
	if (spatial_dimension == 2) extype = "edge2d3";
	else                        extype = "edge3";
	break ;

      case shards::ShellLine<2>::key : extype.assign( "shellline2d2" ); break ;
      case shards::ShellLine<3>::key : extype.assign( "shellline2d3" ); break ;

      case shards::Triangle<3>::key :
	extype = "tri3";
	break ;
      case shards::Triangle<4>::key :
        extype = "tri4";
	break ;
      case shards::Triangle<6>::key :
        extype = "tri6";
	break ;

      case shards::ShellTriangle<3>::key : extype.assign( "trishell3" ); break ;
      case shards::ShellTriangle<6>::key : extype.assign( "trishell6" ); break ;

      case shards::Quadrilateral<4>::key :
        extype = "quad4";
	break ;
      case shards::Quadrilateral<8>::key :
        extype = "quad8";
	break ;
      case shards::Quadrilateral<9>::key :
        extype = "quad9";
	break ;

      case shards::ShellQuadrilateral<4>::key : extype.assign( "shell4" ); break ;
      case shards::ShellQuadrilateral<8>::key : extype.assign( "shell8" ); break ;
      case shards::ShellQuadrilateral<9>::key : extype.assign( "shell9" ); break ;

      case shards::Tetrahedron< 4>::key : extype.assign( "tetra4" ); break ;
      case shards::Tetrahedron<8>::key : extype.assign( "tetra8" ); break ;
      case shards::Tetrahedron<10>::key : extype.assign( "tetra10" ); break ;

      case shards::Pyramid< 5>::key : extype.assign( "pyramid5" ); break ;
      case shards::Pyramid<13>::key : extype.assign( "pyramid13" ); break ;
      case shards::Pyramid<14>::key : extype.assign( "pyramid14" ); break ;

      case shards::Wedge< 6>::key : extype.assign( "wedge6" ); break ;
      case shards::Wedge<15>::key : extype.assign( "wedge15" ); break ;
      case shards::Wedge<18>::key : extype.assign( "wedge18" ); break ;

      case shards::Hexahedron< 8>::key : extype.assign( "hex8" ); break ;
      case shards::Hexahedron<20>::key : extype.assign( "hex20" ); break ;
      case shards::Hexahedron<27>::key : extype.assign( "hex27" ); break ;

      default:
	std::ostringstream oss;
	oss << "stk::io::map_topology_to_ioss( '" << cell_top->name
	    << "' ) ERROR unmapped topology" << std::endl ;
	throw std::runtime_error(oss.str());
      }

      return extype ;
    }

    void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::fem::FEMMetaData &meta,
				  stk::mesh::EntityRank type)
    {
      internal_part_processing(entity, meta.get_meta_data(meta), type);
    }
    void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::fem::FEMMetaData &meta,
				  stk::mesh::EntityRank type)
    {
      internal_part_processing(entity, meta.get_meta_data(meta), type);
    }

    void internal_part_processing(Ioss::GroupingEntity *entity, stk::mesh::MetaData &meta,
				  stk::mesh::EntityRank type)
    {
      if (include_entity(entity)) {
	stk::mesh::Part & part = stk::mesh::declare_part(meta, entity->name(), type);
	stk::io::put_io_part_attribute(part, entity);
      }
    }

    void internal_part_processing(Ioss::EntityBlock *entity, stk::mesh::MetaData &meta,
				  stk::mesh::EntityRank type)
    {
      if (include_entity(entity)) {
        stk::mesh::fem::FEMMetaData * fem_meta = const_cast<stk::mesh::fem::FEMMetaData *>(meta.get_attribute<stk::mesh::fem::FEMMetaData>());
        //const stk::mesh::fem::FEMMetaData * fem_meta = meta.get_attribute<stk::mesh::fem::FEMMetaData>();
        stk::mesh::Part * part = NULL;
        if( fem_meta )
          part = &fem_meta->declare_part(entity->name(), type);
        else
          part = &stk::mesh::declare_part(meta, entity->name(), type);
	stk::io::put_io_part_attribute(*part, entity);

	const Ioss::ElementTopology *topology = entity->topology();
	const CellTopologyData * const cell_topology = map_topology_ioss_to_cell(topology);
	/// \todo IMPLEMENT Determine whether application can work
	/// with this topology type... Perhaps map_topology_ioss_to_cell only
	/// returns a valid topology if the application has registered
        /// that it can handle that specific topology.

        if (cell_topology != NULL) {
          if( fem_meta )
          {
            const stk::mesh::fem::CellTopology cell_topo(cell_topology);
            stk::mesh::fem::set_cell_topology(*fem_meta, *part, cell_topo);
          }
          stk::io::set_cell_topology(*part, cell_topology);
        } else {
	  /// \todo IMPLEMENT handle cell_topolgy mapping error...
	}
	stk::io::define_io_fields(entity, Ioss::Field::ATTRIBUTE, *part, type);
      }
    }

    //----------------------------------------------------------------------
    /** Add all stk::Fields on the entities of the specified part_type
     *  on the specified part of the specified role * to the specified
     *  Ioss::GroupingEntity
     */
    void ioss_add_fields(stk::mesh::Part &part,
			 stk::mesh::EntityRank part_type,
			 Ioss::GroupingEntity *entity,
			 const Ioss::Field::RoleType filter_role,
			 bool add_all)
    {
      stk::mesh::MetaData & meta = mesh::MetaData::get(part);
      stk::mesh::Part &universal = meta.universal_part();

      const std::vector<mesh::FieldBase*> &fields = meta.get_fields();

      std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
      while (I != fields.end()) {
	const stk::mesh::FieldBase *f = *I; ++I;
	if (stk::io::is_valid_part_field(f, part_type, part, universal, filter_role, add_all)) {
	  const stk::mesh::FieldBase::Restriction &res = f->restriction(part_type, part);
	  std::pair<std::string, Ioss::Field::BasicType> field_type;
	  get_io_field_type(f, res.stride[0], &field_type);
	  if (field_type.second != Ioss::Field::INVALID) {
	    int entity_size = entity->get_property("entity_count").get_int();
	    const std::string& name = f->name();
	    entity->field_add(Ioss::Field(name, field_type.second, field_type.first,
					  filter_role, entity_size));
	  }
	}
      }
    }

    /**
     * For the given Ioss::GroupingEntity "entity", find all fields that
     * exist on the input database of type "role" and declare them on
     * the give stk::mesh::Part "part". The "part_type" argument
     * specifies the entity type (node, element, ...) that the field
     * should be declared on for this "part"
     *
     * The "role" will typically be either "ATTRIBUTE" or "TRANSIENT"
     */
    void define_io_fields(Ioss::GroupingEntity *entity,
			  Ioss::Field::RoleType role,
			  stk::mesh::Part &part,
			  stk::mesh::EntityRank part_type)
    {
      stk::mesh::MetaData &meta = mesh::MetaData::get(part);

      bool use_cartesian_for_scalar = false;
      if (role == Ioss::Field::ATTRIBUTE)
	use_cartesian_for_scalar = true;
      
      Ioss::NameList names;
      entity->field_describe(role, &names);

      for (Ioss::NameList::const_iterator I = names.begin(); I != names.end(); ++I) {
	/// \todo IMPLEMENT Need a field selection mechanism and a field naming
	/// (ioss_name -> stk::name)  For now, select all and give the
	/// stk field the same name as the ioss field.

	/// Skip the attribute field that is named "attribute"
	if (*I == "attribute" && names.size() > 1)
	  continue;
	
	/// \todo IMPLEMENT Need to determine whether these are
	/// multi-state fields or constant, or interpolated, or ...
	Ioss::Field io_field = entity->get_field(*I);
	declare_ioss_field(meta, part_type, part, io_field, use_cartesian_for_scalar);
      }
    }

    void get_entity_list(Ioss::GroupingEntity *io_entity,
			 stk::mesh::EntityRank part_type,
			 const stk::mesh::BulkData &bulk,
			 std::vector<stk::mesh::Entity*> &entities)
    {
      std::vector<int> ids ;
      io_entity->get_field_data("ids", ids);

      size_t count = ids.size();
      entities.resize(count);

      for(size_t i=0; i<count; ++i) {
	entities[i] = bulk.get_entity( part_type, ids[i] );
      }
    }

    void field_data_from_ioss(const stk::mesh::FieldBase *field,
			      std::vector<stk::mesh::Entity*> &entities,
			      Ioss::GroupingEntity *io_entity,
			      const std::string &io_fld_name)
    {
      /// \todo REFACTOR Need some additional compatability checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != NULL && io_entity->field_exists(io_fld_name)) {
	Ioss::Field io_field = io_entity->get_field(io_fld_name);
	if (field->type_is<double>()) {
	  internal_field_data_from_ioss(io_field, field, entities, io_entity,
					static_cast<double>(1.0));
	} else if (field->type_is<int>()) {
	  internal_field_data_from_ioss(io_field, field, entities, io_entity,
					static_cast<int>(1));
	}
      }
    }

    void field_data_to_ioss(const stk::mesh::FieldBase *field,
			    std::vector<stk::mesh::Entity*> &entities,
			    Ioss::GroupingEntity *io_entity,
			    const std::string &io_fld_name)
    {
      /// \todo REFACTOR Need some additional compatability checks between
      /// Ioss field and stk::mesh::Field; better error messages...

      if (field != NULL && io_entity->field_exists(io_fld_name)) {
	Ioss::Field io_field = io_entity->get_field(io_fld_name);
	if (field->type_is<double>()) {
	  internal_field_data_to_ioss(io_field, field, entities, io_entity,
				      static_cast<double>(1.0));
	} else if (field->type_is<int>()) {
	  internal_field_data_to_ioss(io_field, field, entities, io_entity,
				      static_cast<int>(1));
	}
      }
    }

    //----------------------------------------------------------------------

    // ----------------------------------------------------------------------
    // Returns true if 'entity' should be a 'part' in the analysis mesh.
    // Returns false if the application is only using a subset of the
    // database entities and this entity is not to be used. The
    // "omitted" property is set by the application during parsing or
    // pre-mesh reading time.
    bool include_entity(Ioss::GroupingEntity *entity)
    {
      assert(entity);

      // Check whether entity has "omitted" property...
      bool omitted = (entity->property_exists("omitted")) &&
	(entity->get_property("omitted").get_int() == 1);

      return !omitted;
    }

    namespace {

      void define_face_edge_block(stk::mesh::Part &part,
				  Ioss::GroupingEntity *ef_set,
				  stk::mesh::EntityRank type,
				  size_t ef_count, int spatial_dimension)
      {
        const stk::mesh::EntityRank siderank = side_rank(mesh::MetaData::get(part));
        const stk::mesh::EntityRank edgerank = edge_rank(mesh::MetaData::get(part));
	ThrowRequire(type == siderank || type == edgerank);

  const CellTopologyData *const ef_topology =   stk::io::get_cell_topology(part) ?
                                                stk::io::get_cell_topology(part) :
                                                stk::mesh::fem::FEMMetaData::get(part).get_cell_topology(part).getCellTopologyData();


	if (ef_topology == NULL) {
          std::ostringstream msg ;
	  msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
	  throw std::runtime_error( msg.str() );
	}

	std::string io_topo = map_topology_cell_to_ioss(ef_topology, spatial_dimension);
	std::string element_topo_name = "unknown";

	// Get faceblock parent element topology quantities...
	// Try to decode from part name...
	std::vector<std::string> tokens;
        stk::util::tokenize(part.name(), "_", tokens);
	if (tokens.size() >= 4) {
	  // Name of form: "name_eltopo_facetopo_id" or
	  //               "name_block_id_facetopo_id"
	  // "name" is typically "surface".
	  const Ioss::ElementTopology *element_topo = Ioss::ElementTopology::factory(tokens[tokens.size()-3], true);
	  if (element_topo != NULL) {
	    element_topo_name = element_topo->name();
	  }
	}

	Ioss::EntityBlock *efb = NULL;
	if (type == siderank) {
	  Ioss::FaceBlock *fb = new Ioss::FaceBlock( ef_set->get_database() ,
						     part.name() ,
						     io_topo, element_topo_name, ef_count);
	  assert(((Ioss::FaceSet*)ef_set)->get_face_block(part.name()) == NULL);
	  ((Ioss::FaceSet*)ef_set)->add(fb);
	  efb = fb;
	} else {
	  Ioss::EdgeBlock *eb = new Ioss::EdgeBlock( ef_set->get_database() , part.name() ,
						     io_topo, element_topo_name, ef_count);
	  assert(((Ioss::EdgeSet*)ef_set)->get_edge_block(part.name()) == NULL);
	  ((Ioss::EdgeSet*)ef_set)->add(eb);
	  efb = eb;
	}

	const mesh::Field<double, mesh::ElementNode> *df = get_distribution_factor_field(part);
	if (df != NULL) {
	  int nodes_per_ef = ef_topology->node_count;
	  std::string storage_type = "Real[";
	  storage_type += Ioss::Utils::to_string(nodes_per_ef);
	  storage_type += "]";
	  efb->field_add(Ioss::Field("distribution_factors", Ioss::Field::REAL, storage_type,
				     Ioss::Field::MESH, ef_count));
	}
      }

      void define_face_edge_blocks(stk::mesh::Part &part,
				   const stk::mesh::BulkData &bulk_data,
				   Ioss::GroupingEntity *ef_set,
				   stk::mesh::EntityRank type,
				   int spatial_dimension)
      {
	mesh::MetaData & meta = mesh::MetaData::get(part);
	assert(type == face_rank(meta) || type == edge_rank(meta));

	const stk::mesh::PartVector &blocks = part.subsets();
	if (blocks.size() > 0) {
	  for (size_t j = 0; j < blocks.size(); j++) {
	    mesh::Part & efb_part = *blocks[j];
            mesh::Selector selector = ( meta.locally_owned_part() | meta.globally_shared_part() ) & efb_part;
	    size_t num_ef = count_selected_entities(selector, bulk_data.buckets(type));

	    define_face_edge_block(efb_part, ef_set, type, num_ef, spatial_dimension);
	  }
	} else {
          mesh::Selector selector = ( meta.locally_owned_part() | meta.globally_shared_part() ) & part;
	  size_t num_ef = count_selected_entities(selector, bulk_data.buckets(type));
	  define_face_edge_block(part, ef_set, type, num_ef, spatial_dimension);
	}
      }

      //----------------------------------------------------------------------
      //----------------------------------------------------------------------
      //----------------------------------------------------------------------

      void define_node_block(stk::mesh::Part &part,
			     const stk::mesh::BulkData &bulk,
			     Ioss::Region &io_region)
      {
	//--------------------------------
	// Set the spatial dimension:
	mesh::MetaData & meta = mesh::MetaData::get(part);

	/// \todo REFACTOR The coordinate field would typically be
	/// stored by the app and wouldn't need to be accessed via
	/// string lookup.  App infrastructure is not shown here, so
	/// lookup by string for the example.
	mesh::Field<double, mesh::Cartesian> *coord_field =
	  meta.get_field<stk::mesh::Field<double, mesh::Cartesian> >(std::string("coordinates"));
	assert(coord_field != NULL);
	const mesh::FieldBase::Restriction &res = coord_field->restriction(node_rank(meta), part);
	
	/** \todo REFACTOR  Need a clear way to query dimensions
	 *                  from the field restriction.
	 */
	const int spatial_dim = res.stride[0] ;
	io_region.property_add( Ioss::Property("spatial_dimension", spatial_dim));

	//--------------------------------
	// Create the special universal node block:

	mesh::Selector selector = meta.locally_owned_part() | meta.globally_shared_part();

	size_t num_nodes = count_selected_entities(selector, bulk.buckets(node_rank(meta)));

	const std::string name("nodeblock_1");

	Ioss::NodeBlock * const nb = new Ioss::NodeBlock(io_region.get_database(),
							   name, num_nodes, spatial_dim);
	io_region.add( nb );
      }


      void define_element_block(stk::mesh::Part &part,
				const stk::mesh::BulkData &bulk,
				Ioss::Region &io_region)
      {

	mesh::MetaData & meta = mesh::MetaData::get(part);
        const stk::mesh::EntityRank elem_rank = element_rank(meta);

      const CellTopologyData * const cell_top = 
        stk::io::get_cell_topology(part) ?
        stk::io::get_cell_topology(part) :
        stk::mesh::fem::FEMMetaData::get(part).get_cell_topology(part).getCellTopologyData();

	if (cell_top == NULL) {
          std::ostringstream msg ;
	  msg << " INTERNAL_ERROR: Part " << part.name() << " returned NULL from get_cell_topology()";
	  throw std::runtime_error( msg.str() );
	}

        mesh::Selector selector = ( meta.locally_owned_part() | meta.globally_shared_part() ) & part;
	const size_t num_elems = count_selected_entities( selector, bulk.buckets(elem_rank));

	int spatial_dim = io_region.get_property("spatial_dimension").get_int();

	// Defer the counting of attributes until after we define the
	// element block so we can count them as we add them as fields to
	// the element block
	int attribute_count = 0;
	Ioss::ElementBlock *eb = new Ioss::ElementBlock(io_region.get_database() ,
							part.name() ,
							map_topology_cell_to_ioss(cell_top, spatial_dim) ,
							num_elems , attribute_count);
	io_region.add(eb);

	// Add the attribute fields.
	ioss_add_fields(part, part_primary_entity_rank(part), eb, Ioss::Field::ATTRIBUTE);
      }

      void define_face_set(stk::mesh::Part &part,
			   const stk::mesh::BulkData &bulk,
			   Ioss::Region &io_region)
      {
        const stk::mesh::EntityRank fa_rank = face_rank(mesh::MetaData::get(part));

	bool create_faceset = true;
	if (part.subsets().empty()) {
	  // Only define a faceset for this part if its superset part is
	  // not a face-containing part..  (i.e., this part is not a subset part
	  // in a surface...)
	  const stk::mesh::PartVector &supersets = part.supersets();
	  for (size_t i=0; i < supersets.size(); i++) {
	    if (is_part_io_part(*supersets[i]) && supersets[i]->primary_entity_rank() == fa_rank) {
	      create_faceset = false;
	      break;
	    }
	  }
	}
	if (create_faceset) {
	  Ioss::FaceSet * const fs = new Ioss::FaceSet(io_region.get_database(), part.name());

	  io_region.add(fs);
	  int spatial_dim = io_region.get_property("spatial_dimension").get_int();
          define_face_edge_blocks(part, bulk, fs, fa_rank, spatial_dim);
	}
      }

      void define_edge_set(stk::mesh::Part &part,
			   const stk::mesh::BulkData &bulk,
			   Ioss::Region &io_region)
      {
        const stk::mesh::EntityRank edge_rank = stk::io::edge_rank(mesh::MetaData::get(part));

	bool create_edgeset = true;
	if (part.subsets().empty()) {
	  // Only define a edgeset for this part if its superset part is
	  // not a edge-containing part..  (i.e., this part is not a subset part
	  // in a surface...)
	  const stk::mesh::PartVector &supersets = part.supersets();
	  for (size_t i=0; i < supersets.size(); i++) {
	    if (is_part_io_part(*supersets[i]) && supersets[i]->primary_entity_rank() == edge_rank) {
	      create_edgeset = false;
	      break;
	    }
	  }
	}
	if (create_edgeset) {
	  Ioss::EdgeSet * const fs = new Ioss::EdgeSet(io_region.get_database(), part.name());

	  io_region.add(fs);
	  int spatial_dim = io_region.get_property("spatial_dimension").get_int();
	  define_face_edge_blocks(part, bulk, fs, edge_rank, spatial_dim);
	}
      }

      void define_node_set(stk::mesh::Part &part,
			   const stk::mesh::BulkData &bulk,
			   Ioss::Region &io_region)
      {
	mesh::MetaData & meta = mesh::MetaData::get(part);

        mesh::Selector selector = ( meta.locally_owned_part() | meta.globally_shared_part() ) & part;

	const size_t num_nodes =
	  count_selected_entities(selector, bulk.buckets(node_rank(meta)));

	Ioss::NodeSet * const ns =
	  new Ioss::NodeSet( io_region.get_database(), part.name(), num_nodes);
	io_region.add(ns);
      }
    } // namespace <blank>

    void define_output_db(Ioss::Region & io_region ,
			  const mesh::BulkData &bulk_data,
			  const Ioss::Region *input_region)
    {
      const mesh::MetaData & meta_data = mesh::MetaData::get(bulk_data);

      const stk::mesh::EntityRank no_rank = node_rank(meta_data);
      const stk::mesh::EntityRank el_rank = element_rank(meta_data);
      const stk::mesh::EntityRank fa_rank = face_rank(meta_data);
      const stk::mesh::EntityRank ed_rank = edge_rank(meta_data);

      io_region.begin_mode( Ioss::STATE_DEFINE_MODEL );

      define_node_block(meta_data.universal_part(), bulk_data, io_region);

      // All parts of the meta data:
      const mesh::PartVector & all_parts = meta_data.get_parts();
      for (mesh::PartVector::const_iterator i = all_parts.begin();
	   i != all_parts.end(); ++i) {

	mesh::Part * const part = *i ;

	if (is_part_io_part(*part)) {
	  if (invalid_rank(part->primary_entity_rank()))
	    continue;
          else if (part->primary_entity_rank() == no_rank)
	    define_node_set(*part, bulk_data, io_region);
          else if (part->primary_entity_rank() == el_rank)
	    define_element_block(*part, bulk_data, io_region);
          else if (part->primary_entity_rank() == fa_rank)
	    define_face_set(*part, bulk_data, io_region);
          else if (part->primary_entity_rank() == ed_rank)
	    define_edge_set(*part, bulk_data, io_region);
	}
      }

      if (input_region != NULL)
	io_region.synchronize_id_and_name(input_region, true);

      io_region.end_mode( Ioss::STATE_DEFINE_MODEL );
    }

    //----------------------------------------------------------------------

    namespace {

      size_t get_entities(stk::mesh::Part &part,
			  const stk::mesh::BulkData &bulk,
			  std::vector<mesh::Entity*> &entities)
      {
	mesh::MetaData & meta = mesh::MetaData::get(part);
	mesh::EntityRank type = part_primary_entity_rank(part);
	if (invalid_rank(type))
	  type = node_rank(meta);

        mesh::Selector selector = part & ( meta.locally_owned_part() | meta.globally_shared_part() );
	get_selected_entities(selector, bulk.buckets(type), entities);
	return entities.size();
      }

      void write_side_data_to_ioss( Ioss::GroupingEntity & io ,
				    mesh::Part * const part ,
				    const mesh::BulkData & bulk_data )
      {
        //std::cout << "tmp write_side_data_to_ioss part= " << part->name() << std::endl;
	const mesh::MetaData & meta_data = mesh::MetaData::get(*part);

	std::vector<mesh::Entity *> sides ;
	size_t num_sides = get_entities(*part, bulk_data, sides);

	std::vector<int> side_ids(num_sides);
	std::vector<int> elem_side_ids(num_sides*2);

	stk::mesh::EntityRank elem_rank = element_rank(meta_data);
	for(size_t i=0; i<num_sides; ++i) {

	  const mesh::Entity &side = *sides[i] ;
	  const mesh::PairIterRelation side_elem = side.relations( elem_rank );

	  // Which element to use?
	  // Any locally owned element that has the "correct" orientation

	  const size_t num_side_elem = side_elem.size();

	  const mesh::Relation *rel = NULL ;

	  for ( size_t j = 0 ; j < num_side_elem && ! rel ; ++j ) {
	    const mesh::Entity & elem = *side_elem[j].entity();

	    if ( elem.bucket().member( meta_data.locally_owned_part() ) &&
		 (num_side_elem == 1 || element_side_polarity(elem, side, side_elem[j].identifier())) ) {
	      rel = &side_elem[j];
	    }
	  }

	  if (rel == NULL) { // no suitable element found
	    std::ostringstream oss;
	    oss << "ERROR, no suitable element found";
	    throw std::runtime_error(oss.str());
	  }

	  side_ids[i] = side.identifier();

	  elem_side_ids[i*2]   = rel->entity()->identifier();
	  elem_side_ids[i*2+1] = rel->identifier() + 1; // Ioss is 1-based, mesh is 0-based.
	}

	const size_t num_ids_written = io.put_field_data("ids", side_ids);
	const size_t num_side_written = io.put_field_data("element_side",elem_side_ids);

	if ( num_sides != num_ids_written || num_sides != num_side_written ) {
	  std::ostringstream msg ;

	  msg << "stk::io::write_side_data_to_ioss FAILED for " ;
	  msg << io.name();
	  msg << " in Ioss::GroupingEntity::put_field_data:" ;
	  msg << " num_sides = " << num_sides ;
	  msg << " , num_ids_written = " << num_ids_written ;
	  msg << " , num_side_written = " << num_side_written ;
	  throw std::runtime_error( msg.str() );
	}

	const mesh::Field<double, mesh::ElementNode> *df = get_distribution_factor_field(*part);
	if (df != NULL) {
	  field_data_to_ioss(df, sides, &io, "distribution_factors");
	}
      }


      //----------------------------------------------------------------------
      void output_node_block(Ioss::NodeBlock &nb,
			     stk::mesh::Part &part,
			     const stk::mesh::BulkData &bulk)
      {
	//----------------------------------
	// Exactly one node block to obtain the nodal coordinates and ids:
	// Note that the "ids" field of the nodes needs to be written
	// before any other bulk data that uses node ids since it sets up
	// the global->local mapping of nodes for the output database.
	// Similarly for the element "ids" field related to bulk data
	// using element ids.
	std::vector<mesh::Entity *> nodes ;
	size_t num_nodes = get_entities(part, bulk, nodes);

	std::vector<int> node_ids(num_nodes);
	for(size_t i=0; i<num_nodes; ++i) {
	  const mesh::Entity & node = * nodes[i] ;
	  node_ids[i] = node.identifier();
	}

	size_t num_ids_written = nb.put_field_data("ids", node_ids);
	if ( num_nodes != num_ids_written) {
	  std::ostringstream msg ;
	  msg << " FAILED in Ioss::NodeBlock::put_field_data:" ;
	  msg << " num_nodes = " << num_nodes ;
	  msg << " , num_ids_written = " << num_ids_written ;
	  throw std::runtime_error( msg.str() );
	}

	/// \todo REFACTOR The coordinate field would typically be
	  /// stored by the app and wouldn't need to be accessed via
	  /// string lookup.  App infrastructure is not shown here, so
	  /// lookup by string for the example.
	  const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
	  mesh::Field<double, mesh::Cartesian> *coord_field =
	    meta_data.get_field<stk::mesh::Field<double, mesh::Cartesian> >(std::string("coordinates"));
	  assert(coord_field != NULL);
	  field_data_to_ioss(coord_field, nodes, &nb, "mesh_model_coordinates");
      }

      void output_element_block(Ioss::ElementBlock *block,
				const stk::mesh::BulkData &bulk)
      {
	const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
	const std::string& name = block->name();
	mesh::Part* part = meta_data.get_part(name);

	assert(part != NULL);
	std::vector<mesh::Entity *> elements;
	size_t num_elems = get_entities(*part, bulk, elements);

	const CellTopologyData * cell_topo = 
              stk::io::get_cell_topology(*part) ?
              stk::io::get_cell_topology(*part) :
              stk::mesh::fem::FEMMetaData::get(*part).get_cell_topology(*part).getCellTopologyData();
	if (cell_topo == NULL) {
          std::ostringstream msg ;
	  msg << " INTERNAL_ERROR: Part " << part->name() << " returned NULL from get_cell_topology()";
	  throw std::runtime_error( msg.str() );
	}
	int nodes_per_elem = cell_topo->node_count;

	std::vector<int> elem_ids(num_elems);
	std::vector<int> connectivity(num_elems*nodes_per_elem);

	stk::mesh::EntityRank no_rank = node_rank(meta_data);
	for (size_t i = 0; i < num_elems; ++i) {

	  elem_ids[i] = elements[i]->identifier();

	  const mesh::PairIterRelation elem_nodes = elements[i]->relations(no_rank);

	  for (int j = 0; j < nodes_per_elem; ++j) {
	    connectivity[i * nodes_per_elem+j] = elem_nodes[j].entity()->identifier();
	  }
	}

	const size_t num_ids_written = block->put_field_data("ids", elem_ids);
	const size_t num_con_written = block->put_field_data("connectivity", connectivity);

	if ( num_elems != num_ids_written || num_elems != num_con_written ) {
	  std::ostringstream msg ;
	  msg << " FAILED in Ioss::ElementBlock::put_field_data:" << std::endl ;
	  msg << "  num_elems = " << num_elems << std::endl ;
	  msg << "  num_ids_written = " << num_ids_written << std::endl ;
	  msg << "  num_connectivity_written = " << num_con_written << std::endl ;
	  throw std::runtime_error( msg.str() );
	}

	stk::mesh::EntityRank elem_rank = element_rank(meta_data);
	const std::vector<mesh::FieldBase *> &fields = meta_data.get_fields();
	std::vector<mesh::FieldBase *>::const_iterator I = fields.begin();
	while (I != fields.end()) {
	  const mesh::FieldBase *f = *I ; ++I ;
	  const Ioss::Field::RoleType *role = stk::io::get_field_role(*f);
	  if (role != NULL && *role == Ioss::Field::ATTRIBUTE) {
	    const mesh::FieldBase::Restriction &res = f->restriction(elem_rank, *part);
	    if (res.stride[0] > 0) {
	      stk::io::field_data_to_ioss(f, elements, block, f->name());
	    }
	  }
	}
      }

      void output_node_set(Ioss::NodeSet *ns, const stk::mesh::BulkData &bulk)
      {
	const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
	const std::string& name = ns->name();
	stk::mesh::Part* part = meta_data.get_part(name);
	assert(part != NULL);

	std::vector<stk::mesh::Entity *> nodes ;
	size_t num_nodes = get_entities(*part, bulk, nodes);

	std::vector<int> node_ids(num_nodes);

	for(size_t i=0; i<num_nodes; ++i) {
	  const stk::mesh::Entity & node = * nodes[i] ;
	  node_ids[i] = node.identifier();
	}
	size_t num_ids_written = ns->put_field_data("ids", node_ids);
	if ( num_nodes != num_ids_written ) {
	  std::ostringstream msg ;
	  msg << " FAILED in Ioss::NodeSet::put_field_data:"
	      << " num_nodes = " << num_nodes
	      << ", num_ids_written = " << num_ids_written;
	  throw std::runtime_error( msg.str() );
	}

	stk::mesh::Field<double> *df_field =
	  meta_data.get_field<stk::mesh::Field<double> >("distribution_factors");
	if (df_field != NULL) {
	  stk::io::field_data_to_ioss(df_field, nodes, ns, "distribution_factors");
	}
      }

      void output_face_set(Ioss::FaceSet *fs,
			   const stk::mesh::BulkData &bulk)
      {
        //std::cout << "tmp output_face_set " << std::endl;
	const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
	size_t block_count = fs->block_count();
	for (size_t i=0; i < block_count; i++) {
	  Ioss::EntityBlock *block = fs->get_block(i);
	  if (stk::io::include_entity(block)) {
	    stk::mesh::Part * const part = meta_data.get_part(block->name());
	    stk::io::write_side_data_to_ioss(*block, part, bulk);
	  }
	}
      }

      void output_edge_set(Ioss::EdgeSet *es,
			   const stk::mesh::BulkData &bulk)
      {
        //std::cout << "tmp output_edge_set " << std::endl;
	const stk::mesh::MetaData & meta_data = mesh::MetaData::get(bulk);
	size_t block_count = es->block_count();
	for (size_t i=0; i < block_count; i++) {
	  Ioss::EntityBlock *block = es->get_block(i);
	  if (stk::io::include_entity(block)) {
	    stk::mesh::Part * const part = meta_data.get_part(block->name());
	    stk::io::write_side_data_to_ioss(*block, part, bulk);
	  }
	}
      }
    } // namespace <blank>

    void write_output_db(Ioss::Region& io_region,
			 const stk::mesh::BulkData& bulk)
    {
	    const stk::mesh::MetaData & meta = mesh::MetaData::get(bulk);

      io_region.begin_mode( Ioss::STATE_MODEL );

      Ioss::NodeBlock & nb = *io_region.get_node_blocks()[0];
      output_node_block(nb, meta.universal_part(), bulk);

      //----------------------------------
      const Ioss::ElementBlockContainer& elem_blocks = io_region.get_element_blocks();
      for(Ioss::ElementBlockContainer::const_iterator it = elem_blocks.begin();
	  it != elem_blocks.end(); ++it) {
	output_element_block(*it, bulk);
      }

      //----------------------------------
      const Ioss::NodeSetContainer& node_sets = io_region.get_nodesets();
      for(Ioss::NodeSetContainer::const_iterator it = node_sets.begin();
	  it != node_sets.end(); ++it) {
	output_node_set(*it, bulk);
      }

      //----------------------------------
      const Ioss::FaceSetContainer& face_sets = io_region.get_facesets();
      for(Ioss::FaceSetContainer::const_iterator it = face_sets.begin();
	  it != face_sets.end(); ++it) {
	output_face_set(*it, bulk);
      }

      //----------------------------------
      const Ioss::EdgeSetContainer& edge_sets = io_region.get_edgesets();
      for(Ioss::EdgeSetContainer::const_iterator it = edge_sets.begin();
	  it != edge_sets.end(); ++it) {
	output_edge_set(*it, bulk);
      }

      io_region.end_mode( Ioss::STATE_MODEL );
    }

    //----------------------------------------------------------------------
    bool is_part_io_part(stk::mesh::Part &part)
    {
      return NULL != part.attribute<Ioss::GroupingEntity>();
    }

    const stk::mesh::Field<double, stk::mesh::ElementNode> *get_distribution_factor_field(const stk::mesh::Part &p)
    {
      return p.attribute<stk::mesh::Field<double, stk::mesh::ElementNode> >();
    }

    void set_distribution_factor_field(stk::mesh::Part &p,
				       const stk::mesh::Field<double, stk::mesh::ElementNode> &df_field)
    {
      stk::mesh::MetaData &m = mesh::MetaData::get(p);
      m.declare_attribute_no_delete(p,&df_field);
    }

    const Ioss::Field::RoleType* get_field_role(const stk::mesh::FieldBase &f)
    {
      return f.attribute<Ioss::Field::RoleType>();
    }

    void set_field_role(stk::mesh::FieldBase &f, const Ioss::Field::RoleType &role)
    {
      Ioss::Field::RoleType *my_role = new Ioss::Field::RoleType(role);
      stk::mesh::MetaData &m = mesh::MetaData::get(f);
      const Ioss::Field::RoleType *check = m.declare_attribute_with_delete(f, my_role);
      if ( check != my_role ) {
        if (*check != *my_role) {
          std::ostringstream msg ;
	  msg << " FAILED in IossBridge -- set_field_role:"
	      << " The role type had already been set to " << *check
	      << ", so it is not possible to change it to " << *my_role;
	  throw std::runtime_error( msg.str() );
        }
        delete my_role;
      }
    }

  }//namespace io
}//namespace stk

