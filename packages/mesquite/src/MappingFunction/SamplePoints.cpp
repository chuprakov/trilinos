/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2006 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

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


/** \file SamplePoints.cpp
 *  \brief 
 *  \author Jason Kraftcheck 
 */

#include "Mesquite.hpp"
#include "SamplePoints.hpp"
#include "TopologyInfo.hpp"
#include <assert.h>
#include <stdlib.h>

namespace Mesquite {

const unsigned char CORNERS = 1;
const unsigned char MIDEDGE = (CORNERS << 1);
const unsigned char MIDFACE = (CORNERS << 2);
const unsigned char MIDELEM = (CORNERS << 3);

SamplePoints::SamplePoints( bool do_default )
{
  clear_all();
  if (do_default) {
    sampleBits[POLYGON      ] = CORNERS;
    sampleBits[TRIANGLE     ] = MIDFACE;
    sampleBits[QUADRILATERAL] = CORNERS;
    sampleBits[POLYHEDRON   ] = 0;
    sampleBits[TETRAHEDRON  ] = MIDELEM;
    sampleBits[HEXAHEDRON   ] = CORNERS;
    sampleBits[PRISM        ] = CORNERS;
    sampleBits[PYRAMID      ] = CORNERS;
    sampleBits[SEPTAHEDRON  ] = CORNERS;
  }
}

void SamplePoints::clear_all()
  { for (size_t i = 0; i < MIXED; sampleBits[i++] = 0); }

SamplePoints::SamplePoints( bool corners,
                            bool midedge,
                            bool midface,
                            bool midelem)
{
  size_t i;
  unsigned char bits = 0;
  if (corners) bits |= CORNERS;
  if (midedge) bits |= MIDEDGE;
  if (midface) bits |= MIDFACE;
  for (i = 0; i < TETRAHEDRON; ++i)
    sampleBits[i] = bits;
  if (midelem) bits |= MIDELEM;
  for (i = TETRAHEDRON; i < MIXED; ++i)
    sampleBits[i] = bits;
}
  
unsigned SamplePoints::num_sample_points( EntityTopology type, unsigned bits )
{
  const unsigned dim = TopologyInfo::dimension( type );
  unsigned count = !!(bits & (1 << dim));
  for (unsigned i = 0; i < dim; ++i)
    if (bits & (1 << i))
      count += TopologyInfo::adjacent( type, i );
    // special case for pyramids: we never want to evaluate at
    // the apex of the pyramid
  if (type == PYRAMID && bits&1)
    --count;
  return count;
}

void SamplePoints::location_from_sample_number( EntityTopology element_type,
                                           unsigned sample_topology_bits,
                                           unsigned sample_point_number,
                                           unsigned& entity_dimension,
                                           unsigned& entity_number )
{
  const unsigned dim = TopologyInfo::dimension( element_type );
  for (unsigned i = 0; i < dim; ++i) {
    if ( sample_topology_bits & (1 << i)) {
      unsigned n =  TopologyInfo::adjacent( element_type, i );
        // Skip apex for pyramid elements
      if (element_type == PYRAMID && i == 0)
        --n;
      if (n > sample_point_number) {
        entity_number = sample_point_number;
        entity_dimension = i;
        return;
      }
      sample_point_number -= n;
    }
  }
  
  assert( sample_point_number == 0 );
  entity_dimension = dim;
  entity_number = 0;
}

unsigned SamplePoints::sample_number_from_location( EntityTopology element_type,
                                           unsigned sample_topology_bits,
                                           unsigned entity_dimension,
                                           unsigned entity_number )
{
  const unsigned dim = TopologyInfo::dimension( element_type );
  assert( entity_dimension <= dim );
  if (!(sample_topology_bits & (1 << entity_dimension)))
    return (unsigned)-1; 
  
  unsigned result = (entity_dimension == dim) ? 0 : entity_number;
  for (unsigned i = 0; i < entity_dimension; ++i) 
    if (sample_topology_bits & (1 << i)) 
      result += TopologyInfo::adjacent( element_type, i );
  
    // skip apex for pyramids
  if (element_type == PYRAMID && entity_dimension > 0)
    --result;
    
  return result;
}

} // namespace Mesquite
