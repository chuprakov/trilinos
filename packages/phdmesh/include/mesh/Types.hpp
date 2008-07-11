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

#include <util/Basics.hpp>
#include <util/Span.hpp>
#include <util/Dimension.hpp>

namespace phdmesh {

//----------------------------------------------------------------------
// Mesh meta-data

class Schema ;  // Meta-data description of a mesh
class Part ;    // Defined subset of the mesh

template< typename Scalar = void , class Tag1 = DimensionTag ,
                                   class Tag2 = DimensionTag ,
                                   class Tag3 = DimensionTag ,
                                   class Tag4 = DimensionTag ,
                                   class Tag5 = DimensionTag ,
                                   class Tag6 = DimensionTag ,
                                   class Tag7 = DimensionTag >
class Field ;

typedef
Field< void , DimensionTag , DimensionTag ,
              DimensionTag , DimensionTag ,
              DimensionTag , DimensionTag ,
              DimensionTag > FieldBase ;

enum { MaximumFieldDimension = 7 };

//----------------------------------------------------------------------
// Mesh bulk-data

class Mesh ;     // Bulk-data of a mesh
class Kernel ;   // Homogeneous collection of mesh entitities their field data
class Entity ;   // Individual entity within the mesh
class Relation ; // Relation pair of local mesh entities

//----------------------------------------------------------------------
// Supporting types for entity and relation attributes

typedef uint32_type   entity_id_type ;     // Entity identifier type
typedef uint64_type   entity_key_type ;    // Entity key type
typedef uint_ptr_type relation_attr_type ; // Entity relation attribute

/** Extensible definition of types of entities.
 *  The first four types are required to have values
 *  corresponding to the topological rank of the entity type.
 *  The types values must be sequencial and 'EntityTypeEnd'
 *  must be the last value.
 */
enum EntityType {
  Node       = 0 ,
  Edge       = 1 ,
  Face       = 2 ,
  Element    = 3 ,
  Particle   = 4 ,
  Constraint = 5 ,
  EntityTypeEnd = 6 };

/** Number of binary digits used to hold an entity type value */
enum { entity_key_type_digits = 4 };

const char * entity_type_name( EntityType );

//----------------------------------------------------------------------
/** A relation stencil maps entity relationships to ordinals.
 *  If the given relationship is not in the domain of the stencil
 *  then a negative value must be returned.
 */
typedef int ( * relation_stencil_ptr )(
  EntityType from_type ,
  EntityType to_type ,
  unsigned   identifier ,
  unsigned   kind );

//----------------------------------------------------------------------
// Parallel bulk-data

typedef std::pair<Entity*,unsigned> EntityProc ; // Entity-processor pair

typedef std::vector< EntityProc > EntityProcSet ;

typedef Span< EntityProcSet::const_iterator > EntityProcSpan ;

//----------------------------------------------------------------------

} // namespace phdmesh

#endif

