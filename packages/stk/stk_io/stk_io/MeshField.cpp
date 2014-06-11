/*--------------------------------------------------------------------*/
/*    Copyright 2004, 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#include <stk_io/MeshField.hpp>
#include <stk_io/DbStepTimeInterval.hpp>
#include <stk_io/IossBridge.hpp>        
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase, etc
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field, etc
#include <stk_util/environment/ReportHandler.hpp>

#include <Ioss_GroupingEntity.h>        // for GroupingEntity
#include <Ioss_VariableType.h>          // for VariableType

#include <math.h>
#include <assert.h>

namespace stk {
  namespace io {

MeshField::MeshField(stk::mesh::FieldBase *field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(field),
    m_dbName(db_name),
    m_timeToRead(0.0),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_singleState(true),
    m_isActive(false)
{
  if (db_name == "") {
    m_dbName = field->name();
  }

  ThrowErrorMsgIf(m_timeMatch == LINEAR_INTERPOLATION && m_field->type_is<int>(),
		  "ERROR: Input interpolation field '" << m_field->name()
		  << "' is an integer field.  Only double fields can be interpolated.");
}
    
MeshField::MeshField(stk::mesh::FieldBase &field,
		     const std::string &db_name,
		     TimeMatchOption tmo)
  : m_field(&field),
    m_dbName(db_name),
    m_timeToRead(0.0),
    m_timeMatch(tmo),
    m_oneTimeOnly(false),
    m_isActive(false)
{
  if (db_name == "") {
    m_dbName = field.name();
  }

  ThrowErrorMsgIf(m_timeMatch == LINEAR_INTERPOLATION && m_field->type_is<int>(),
		  "ERROR: Input interpolation field '" << m_field->name()
		  << "' is an integer field.  Only double fields can be interpolated.");
}

MeshField::~MeshField()
{}

MeshField& MeshField::set_read_time(double time_to_read)
{
  m_timeToRead = time_to_read;
  m_timeMatch = SPECIFIED;
  m_oneTimeOnly = true;
  return *this;
}

MeshField& MeshField::set_active()
{
  m_isActive = true;
  return *this;
}
MeshField& MeshField::set_inactive()
{
  m_isActive = false;
  return *this;
}

MeshField& MeshField::set_single_state(bool yesno)
{
  m_singleState = yesno;
  return *this;
}

MeshField& MeshField::set_read_once(bool yesno)
{
  m_oneTimeOnly = yesno;
  return *this;
}

void MeshField::add_part(const stk::mesh::EntityRank rank,
			 Ioss::GroupingEntity *io_entity)
{
  m_fieldParts.push_back(MeshFieldPart(rank,io_entity, m_dbName));
}

bool MeshField::operator==(const MeshField &other) const
{
  return m_field       == other.m_field &&
         m_dbName      == other.m_dbName &&
         m_subsetParts == other.m_subsetParts;
}

void MeshField::restore_field_data(stk::mesh::BulkData &bulk,
				   const stk::io::DBStepTimeInterval &sti)
{
  if (!is_active())
    return;
  
  if (m_timeMatch == CLOSEST || m_timeMatch == SPECIFIED) {
    int step = 0;
    if (m_timeMatch == CLOSEST) {
      step = sti.get_closest_step();
    }
    else if (m_timeMatch == SPECIFIED) {
      DBStepTimeInterval sti2(sti.region, m_timeToRead);
      step = sti2.get_closest_step();
    }
    STKIORequire(step > 0);
    
    sti.region->begin_state(step);

    std::vector<stk::io::MeshFieldPart>::iterator I = m_fieldParts.begin();
    while (I != m_fieldParts.end()) {
      const stk::mesh::EntityRank part_rank = (*I).get_entity_rank();
      Ioss::GroupingEntity *io_entity = (*I).get_io_entity();
      std::vector<stk::mesh::Entity> entity_list;
      stk::io::get_entity_list(io_entity, part_rank, bulk, entity_list);
      
      size_t state_count = m_field->number_of_states();
      stk::mesh::FieldState state = m_field->state();
      // If the multi-state field is not "set" at the newest state, then the user has
      // registered the field at a specific state and only that state should be input.
      if(m_singleState || state_count == 1 || state != stk::mesh::StateNew) {
	stk::io::field_data_from_ioss(bulk, m_field, entity_list, io_entity, m_dbName);
      } else {
	stk::io::multistate_field_data_from_ioss(bulk, m_field, entity_list, io_entity, m_dbName, state_count);
      }

      if (m_oneTimeOnly) {
	(*I).release_field_data();
      }
      ++I;
    }
    sti.region->end_state(step);
  }
  else if (m_timeMatch == LINEAR_INTERPOLATION) {
    std::vector<stk::io::MeshFieldPart>::iterator I = m_fieldParts.begin();
    while (I != m_fieldParts.end()) {

      // Get data at beginning of interval...
      std::vector<double> values;
      (*I).get_interpolated_field_data(sti, values);
      
      size_t state_count = m_field->number_of_states();
      stk::mesh::FieldState state = m_field->state();

      // Interpolation only handles single-state fields currently.
      STKIORequire(m_singleState || state_count == 1 || state != stk::mesh::StateNew);

      Ioss::GroupingEntity *io_entity = (*I).get_io_entity();
      const Ioss::Field &io_field = io_entity->get_fieldref(m_dbName);
      size_t field_component_count = io_field.transformed_storage()->component_count();

      std::vector<stk::mesh::Entity> entity_list;
      const stk::mesh::EntityRank part_rank = (*I).get_entity_rank();
      stk::io::get_entity_list(io_entity, part_rank, bulk, entity_list);
      
      for (size_t i=0; i < entity_list.size(); ++i) {
	if (bulk.is_valid(entity_list[i])) {
	  double *fld_data = static_cast<double*>(stk::mesh::field_data(*m_field, entity_list[i]));
	  if (fld_data !=NULL) {
	    for(size_t j=0; j<field_component_count; ++j) {
	      fld_data[j] = values[i*field_component_count+j];
	    }
	  }
	}
      }
      if (m_oneTimeOnly) {
	(*I).release_field_data();
      }
      ++I;
    }
  }
  if (m_oneTimeOnly) {
    set_inactive();
  }
}

void MeshFieldPart::release_field_data()
{
  m_preStep = 0;
  m_postStep = 0;
  std::vector<double>().swap(m_preData);
  std::vector<double>().swap(m_postData);
}

void MeshFieldPart::load_field_data(const DBStepTimeInterval &sti)
{
  // Use cached data if possible; avoid reading from disk...
  
  if (sti.exists_before && m_preStep != sti.s_before) {
    assert(sti.s_before > 0);

    if (sti.s_before == m_postStep) {
      m_preData.swap(m_postData);
      std::swap(m_preStep, m_postStep);
    }
    else {
      // See if postStep can use my current data...
      if (sti.exists_after && sti.s_after == m_preStep) {
	m_postData.swap(m_preData);
	std::swap(m_postStep, m_preStep);
      }
      m_preStep = sti.s_before;
      sti.region->begin_state(m_preStep);
      m_ioEntity->get_field_data(m_dbName, m_preData);
      sti.region->end_state(m_preStep);
    }
  }

  if (sti.exists_after && m_postStep != sti.s_after) {
    m_postStep = sti.s_after;
    assert(m_postStep > 0);

    if (m_preStep == m_postStep) {
      m_postData.resize(0);
      m_postData.reserve(m_preData.size());
      std::copy(m_preData.begin(), m_preData.end(), std::back_inserter(m_postData));
    }
    else {
      sti.region->begin_state(m_postStep);
      m_ioEntity->get_field_data(m_dbName, m_postData);
      sti.region->end_state(m_postStep);
    }
  }
}

void MeshFieldPart::get_interpolated_field_data(const DBStepTimeInterval &sti, std::vector<double> &values)
{
  load_field_data(sti);
  size_t values_size = sti.exists_before ? m_preData.size() : m_postData.size();
  values.resize(0);
  values.reserve(values_size);

  if (sti.exists_before && !sti.exists_after) {
    std::copy(m_preData.begin(), m_preData.end(), std::back_inserter(values));
  }
  else if (!sti.exists_before && sti.exists_after) {
    std::copy(m_postData.begin(), m_postData.end(), std::back_inserter(values));
  }
  else {
    assert(sti.exists_before && sti.exists_after);
    if (sti.s_after == sti.s_before) {
      // No interpolation. preData and postData contain the same step.
      std::copy(m_preData.begin(), m_preData.end(), std::back_inserter(values));
    }
    else {
      // Interpolate
      double tb = sti.region->get_state_time(sti.s_before);
      double ta = sti.region->get_state_time(sti.s_after);
      double delta =  ta - tb;
      double frac  = (sti.t_analysis - tb) / delta;

      assert(m_preData.size() == m_postData.size());
      for (size_t i=0; i < m_preData.size(); i++) {
	values.push_back((1.0 - frac) * m_preData[i] + frac * m_postData[i]);
      }
    }
  }
}

}}  // Close namespaces
