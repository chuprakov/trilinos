/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Sandia National Laboratories.  Developed at the
    University of Wisconsin--Madison under SNL contract number
    624796.  The U.S. Government and the University of Wisconsin
    retain certain rights to this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    (2006) kraftche@cae.wisc.edu
   
  ***************************************************************** */


/** \file DomainSurfaceOrientation.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "DomainSurfaceOrientation.hpp"
#include "SamplePoints.hpp"
#include "MsqMatrix.hpp"
#include "MsqError.hpp"
#include "Vector3D.hpp"
#include "PatchData.hpp"
#include "MappingFunction.hpp"
#include "WeightCalculator.hpp"
#include "ElementQM.hpp"

namespace Mesquite {
int DomainSurfaceOrientation::get_negate_flag( ) const { return 1; }

msq_std::string DomainSurfaceOrientation::get_name() const
  { return msq_std::string("DomainSurfaceOrientation"); }

void DomainSurfaceOrientation::get_evaluations( PatchData& pd,
                                      msq_std::vector<size_t>& handles,
                                      bool free,
                                      MsqError& err )
{
  handles.clear();
  msq_std::vector<size_t> elems;
  ElementQM::get_element_evaluations( pd, elems, free, err ); MSQ_ERRRTN(err);
  for (msq_std::vector<size_t>::iterator i = elems.begin(); i != elems.end(); ++i)
  {
    EntityTopology type = pd.element_by_index( *i ).get_element_type();
    if (TopologyInfo::dimension(type) == 2) {
      unsigned num_samples = samplePts->num_sample_points( type );
      for (unsigned j = 0; j < num_samples; ++j)
        handles.push_back( handle(j, *i) );
    }
  }
}

void DomainSurfaceOrientation::get_element_evaluations( PatchData& pd,
                                              size_t elem,
                                              msq_std::vector<size_t>& handles,
                                              MsqError& err )
{
  EntityTopology type = pd.element_by_index( elem ).get_element_type();
  if (TopologyInfo::dimension(type) == 3) {
    handles.clear();
    return;
  }
  unsigned num_samples = samplePts->num_sample_points( type );
  handles.resize( num_samples );
  for (unsigned j = 0; j < num_samples; ++j)
    handles[j] = handle(j, elem);
}

bool DomainSurfaceOrientation::evaluate( PatchData& pd, size_t handle, double& value, MsqError& err )
{
  mIndices.clear();
  return evaluate_with_indices( pd, handle, value, mIndices, err );
}

bool DomainSurfaceOrientation::evaluate_with_indices( PatchData& pd,
                                            size_t handle,
                                            double& value,
                                            msq_std::vector<size_t>& indices,
                                            MsqError& err )
{
  unsigned s = ElemSampleQM::sample( handle );
  size_t   e = ElemSampleQM::  elem( handle );
  MsqMeshEntity& elem = pd.element_by_index( e );
  EntityTopology type = elem.get_element_type();
  if (TopologyInfo::dimension( type ) != 2) {
    MSQ_SETERR(err)("Surface Orientation metric is valid only for surface elements", 
                    MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }
  
  const MappingFunction* func = pd.get_mapping_function( type );
  if (!func) {
    MSQ_SETERR(err)( "No mapping function for element type", MsqError::UNSUPPORTED_ELEMENT );
    return false;
  }
  
  unsigned dim, num;
  samplePts->location_from_sample_number( type, s, dim, num );
  const unsigned bits = pd.higher_order_node_bits( e );
  
  Vector3D n;
  indices.clear();
  mDerivs.clear();
  switch (dim) {
    case 0:
      func->derivatives_at_corner( num, bits, indices, mDerivs, err );MSQ_ERRZERO(err);
      pd.get_domain_normal_at_corner( e, num, n, err ); MSQ_ERRZERO(err);
      break;
    case 1:
      func->derivatives_at_mid_edge( num, bits, indices, mDerivs, err ); MSQ_ERRZERO(err);
      pd.get_domain_normal_at_mid_edge( e, num, n, err ); MSQ_ERRZERO(err);
      break;
    case 2:
      func->derivatives_at_mid_elem( bits, indices, mDerivs, err ); MSQ_ERRZERO(err);
      pd.get_domain_normal_at_element( e, n, err ); MSQ_ERRZERO(err);
      break;
    default:
      MSQ_SETERR(err)( MsqError::INTERNAL_ERROR );
      return false;
  }
  
    // Convert from indices into element connectivity list to
    // indices into vertex array in patch data.
  const size_t* conn = elem.get_vertex_index_array();
  for (msq_std::vector<size_t>::iterator i = indices.begin(); i != indices.end(); ++i)
    *i = conn[*i];
  
  std::vector<double>::const_iterator d = mDerivs.begin();
  Vector3D c[2] = { Vector3D(0,0,0), Vector3D(0,0,0) };
  for (size_t i = 0; i < indices.size(); ++i) {
    Vector3D coords = pd.vertex_by_index( indices[i] );
    c[0] += *d * coords; ++d;
    c[1] += *d * coords; ++d;
  }
  
  n.normalize();
  const double m[2] = { c[0] % n, c[1] % n };
  value = m[0]*m[0] + m[1]*m[1];
  
    // apply target weight to value
  const double ck = weightCalc->get_weight( pd, e, samplePts, s, err ); MSQ_ERRZERO(err);
  value *= ck;
  return true;
}

} // namespace Mesquite
