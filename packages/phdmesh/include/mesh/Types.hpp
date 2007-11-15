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

#ifndef phdmesh_Types_hpp
#define phdmesh_Types_hpp

//----------------------------------------------------------------------

#include <limits>
#include <utility>
#include <vector>

#include <util/SpanIter.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

class Schema ;  // Parts and fields of a mesh
class Part ;    // Designated subset or part of the mesh

template< typename T , unsigned NDim = 0 > class Field ;

enum { MaximumFieldDimension = 8 };

//----------------------------------------------------------------------

class Mesh ;    // Kernels and entities of a mesh
class Kernel ;  // Homogeneous collection of mesh entitities
class Entity ;  // Individual entity within the mesh
class Connect ; // Connect pair of local mesh entities

typedef std::pair<Entity*,unsigned> EntityProc ; // Entity-processor pair

typedef std::vector< EntityProc > EntityProcSet ;

typedef SpanIter< EntityProcSet::const_iterator > EntityProcSpan ;

//----------------------------------------------------------------------
/** Types of mesh entities.  Extensible via update to the enumeration */

enum EntityType {
  Node    = 0 ,
  Edge    = 1 ,
  Face    = 2 ,
  Element = 3 ,
  Other   = 4 ,
  EntityTypeMaximum = 5 ,
  EntityTypeMask = 0x0f };

enum { EntityTypeDigits = 4 ,
       EntityIdentifierDigits =
         std::numeric_limits<unsigned long>::digits - EntityTypeDigits };

/** Query text name for entity type */

const char * entity_type_name( EntityType );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

