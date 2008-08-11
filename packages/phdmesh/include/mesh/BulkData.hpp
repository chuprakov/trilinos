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

#ifndef phdmesh_Mesh_hpp
#define phdmesh_Mesh_hpp

//----------------------------------------------------------------------

#include <util/Parallel.hpp>
#include <mesh/Types.hpp>
#include <mesh/Field.hpp>
#include <mesh/Entity.hpp>

//----------------------------------------------------------------------

namespace phdmesh {

//----------------------------------------------------------------------
/** Parallel Heterogeneous Dynamic Mesh.
 *  An dynamic unstructured mesh of mesh entities with
 *  subsets of parts partitioned into homogeneous kernels.
 */

class MeshBulkData {
public:

  ~MeshBulkData();

  /** Construct mesh for the given MeshMetaData, parallel machine, and
   *  with the specified maximum number of entities per kernel.
   */
  MeshBulkData( const MeshMetaData & mesh_meta_data ,
                ParallelMachine parallel , unsigned = 1000);

  const MeshMetaData & mesh_meta_data() const { return m_mesh_meta_data ; }

  ParallelMachine parallel() const { return m_parallel_machine ; }
  unsigned parallel_size()   const { return m_parallel_size ; }
  unsigned parallel_rank()   const { return m_parallel_rank ; }

  //------------------------------------
  
  /** Rotation of states:
   *    StateNM1 <- StateNew
   *    StateNM2 <- StateNM1
   *    StateNM3 <- StateNM2
   *    StateNM3 <- StateNM2
   *  etc.
   */
  void update_state();

  //------------------------------------
  /** All kernels of a given entity type */
  const KernelSet & kernels( unsigned type ) const ;

  /** All entitities of a given entity type */
  const EntitySet & entities( unsigned type ) const ;

  Entity * get_entity( entity_key_type ,
                       const char * required_by = NULL ) const ;

  //------------------------------------
  /** Create or retrieve entity of the given type and id.
   *  If the entity is created it is assumed to be locally owned.
   */
  Entity & declare_entity( entity_key_type ,
                           const std::vector<Part*> & ,
                           int parallel_owner = -1 /* default to local */ );

  /** Change the entity's part membership by adding and/or removing parts */
  void change_entity_parts( Entity & ,
                            const std::vector<Part*> & add_parts ,
                            const std::vector<Part*> & remove_parts = 
                                  std::vector<Part*>() );

  void change_entity_owner( Entity & , unsigned );

  void change_entity_identifier( Entity & , entity_id_type );

  /** Declare a relation and its converse between entities in the same mesh.
   *  This mapping ( e_from , identifier , kind ) -> e_to  must be unique.
   */
  void declare_relation( Entity & e_from ,
                         Entity & e_to ,
                         const unsigned identifier ,
                         const unsigned kind = 0 );

  void declare_relation( Entity & , const std::vector<Relation> & );

  /** Remove all relations between two entities. */
  void destroy_relation( Entity & , Entity & , unsigned kind = 0 );

  /** Destroy an entity */
  void destroy_entity( Entity * );

  //------------------------------------
  /** Symmetric parallel relations for shared mesh entities.  */
  const std::vector<EntityProc> & shares() const { return m_shares_all ; }

  /** Asymmetric parallel relations for owner-to-ghosted mesh entities.
   *  Both the source and the destination are fully ordered.
   */
  const std::vector<EntityProc> & ghost_source() const
    { return m_aura_domain ; }

  const std::vector<EntityProc> & ghost_destination() const
    { return m_aura_range ; }

  /** The input must be symmetric, fully ordered, and
   *  each mesh entity must be a member of the 'uses_part'.
   */
  void set_shares( const std::vector<EntityProc> & );

  /** The domain and range inputs must be asymmetric, fully ordered, and
   *  each mesh entity must not be a member of the 'uses_part'.
   *  Domain entities must be owns and range entries must match
   *  the 'owner_field' value.
   */
  void set_ghosting( const std::vector<EntityProc> & ,
                     const std::vector<EntityProc> & );

private:

  MeshBulkData();
  MeshBulkData( const MeshBulkData & );
  MeshBulkData & operator = ( const MeshBulkData & );

  const MeshMetaData      & m_mesh_meta_data ;
  ParallelMachine           m_parallel_machine ;
  unsigned                  m_parallel_size ;
  unsigned                  m_parallel_rank ;
  unsigned                  m_kernel_capacity ;
  KernelSet                 m_kernels[  EntityTypeEnd ];
  EntitySet                 m_entities[ EntityTypeEnd ];
  std::vector<EntityProc>   m_shares_all ;
  std::vector<EntityProc>   m_aura_domain ;
  std::vector<EntityProc>   m_aura_range ;

  /*  Entity modification consequences:
   *  1) Change entity relation => update via part relation => change parts
   *  2) Change parts => update forward relations via part relation
   *                  => update via field relation
   */

  void remove_entity( KernelSet::iterator , unsigned );

  KernelSet::iterator declare_kernel( const EntityType ,
                                      const unsigned ,
                                      const unsigned[] );

  void destroy_kernel( KernelSet::iterator );

  void internal_change_entity_parts( Entity & ,
                                     const PartSet & add_parts ,
                                     const PartSet & remove_parts );

  void internal_propagate_part_changes( Entity & , const PartSet & removed );

  void internal_propagate_relocation( Entity & );
};

//----------------------------------------------------------------------

/** Count entities of each type that are owned by the mesh and
   members of the specified part.

  \param mesh
  \param part
  \param count is an array of length number-of-entity-types (EntityTypeEnd).
*/
void count_entities(
  MeshBulkData & mesh ,
  Part & part ,
  unsigned * const count /* [ EntityTypeEnd ] */ );

/** Count entities of each type that are owned by the mesh and
   members of the specified part.

  \param mesh
  \param part
  \param count is an array of length number-of-entity-types (EntityTypeEnd).
*/
void count_entities(
  MeshBulkData & mesh ,
  const PartSet & parts ,
  unsigned * const count /* [ EntityTypeEnd ] */ );

/** Get all kernels within the given part.
 *  Every kernel will have the part in its superset.
 */
void get_kernels( const KernelSet & , Part & , std::vector<const Kernel*> & );

/** Get all kernels within all of the given parts.
 *  The input PartSet must be properly ordered,
 *  e.g. via the 'phdmesh::order( PartSet & )' function.
 *  Every kernel will have all of the parts in its superset.
 *  It is more efficient to pre-define an intersection part
 *  in the mesh_meta_data and then use the get_kernels function.
 */
void get_kernels_intersect( const KernelSet & ,
                            const PartSet & ,
                            std::vector<const Kernel*> & );

/** Get all kernels within at least one of the given parts.
 *  Every kernel will have at least one of the parts in its superset.
 *  It is more efficient to pre-define a superset part
 *  in the mesh_meta_data and then use the get_kernels function.
 */
void get_kernels_union( const KernelSet & ,
                        const PartSet & ,
                        std::vector<const Kernel*> & );

} // namespace phdmesh

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif

