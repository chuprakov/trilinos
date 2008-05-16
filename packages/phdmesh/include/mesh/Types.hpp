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

class Schema ;  // Parts and fields of a mesh
class Part ;    // Designated subset or part of the mesh

template< typename Scalar = void ,
          class Trait1 = DimensionTraits ,
          class Trait2 = DimensionTraits ,
          class Trait3 = DimensionTraits ,
          class Trait4 = DimensionTraits ,
          class Trait5 = DimensionTraits ,
          class Trait6 = DimensionTraits ,
          class Trait7 = DimensionTraits >
class Field ;

typedef
Field< void , DimensionTraits , DimensionTraits ,
              DimensionTraits , DimensionTraits ,
              DimensionTraits , DimensionTraits ,
              DimensionTraits > FieldBase ;

enum { MaximumFieldDimension = 7 };

//----------------------------------------------------------------------

class Mesh ;    // Kernels and entities of a mesh
class Kernel ;  // Homogeneous collection of mesh entitities
class Entity ;  // Individual entity within the mesh
class Relation ; // Relation pair of local mesh entities

typedef std::pair<Entity*,unsigned> EntityProc ; // Entity-processor pair

typedef std::vector< EntityProc > EntityProcSet ;

typedef Span< EntityProcSet::const_iterator > EntityProcSpan ;

//----------------------------------------------------------------------

typedef uint64_type entity_key_type ;
typedef uint32_type entity_id_type ;

enum {
  entity_key_digits     = std::numeric_limits<entity_key_type>::digits ,
  entity_rank_digits    = 4 ,
  entity_id_digits_max  = entity_key_digits - entity_rank_digits ,
  entity_id_digits_want = std::numeric_limits<entity_id_type>::digits ,
  entity_id_digits      = entity_id_digits_max < entity_id_digits_want ?
                          entity_id_digits_max : entity_id_digits_want
};


enum { end_entity_rank = ((unsigned)        1) << entity_rank_digits };
enum { end_entity_id   = ((entity_key_type) 1) << entity_id_digits };

inline
entity_key_type entity_key( unsigned rank , entity_id_type id )
{
  enum { mask = ( ~((entity_key_type) 0) ) >>
                ( entity_key_digits - entity_id_digits ) };

  return ( ((entity_key_type) rank) << entity_id_digits_max ) | ( id & mask );
}

inline
unsigned entity_rank( entity_key_type key )
{ return key >> entity_id_digits_max ; }

inline
entity_id_type entity_id( entity_key_type key )
{
  enum { mask = ( ~((entity_key_type) 0) ) >>
                ( entity_key_digits - entity_id_digits ) };
  return key & mask ;
}

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

