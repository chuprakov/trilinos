/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_Entity_hpp
#define stk_mesh_Entity_hpp

#include <utility>
#include <vector>

#include <stk_util/util/NamedPair.hpp>
#include <stk_util/util/PairIter.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Relation.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 * \{
 */

//----------------------------------------------------------------------
/** \brief  Span of a sorted relations for a given domain entity.
 *
 *  The span is sorted by
 *  -# range entity rank,
 *  -# relation identifier, and
 *  -# range entity global identifier.
 */
typedef PairIter< std::vector<Relation>::const_iterator > PairIterRelation ;

//----------------------------------------------------------------------

NAMED_PAIR( EntityCommInfo , unsigned , ghost_id , unsigned , proc )

/** \brief  Span of ( communication-subset-ordinal , process-rank ) pairs
 *          for the communication of an entity.
 */
typedef PairIter< std::vector< EntityCommInfo >::const_iterator >
  PairIterEntityComm ;

//----------------------------------------------------------------------
/** \brief  A fundamental unit within the discretization of a problem domain,
 *          including but not limited to nodes, edges, sides, and elements.
 *
 *  Entities are distributed among parallel processors.
 *  A given entity may reside on more than one processor;
 *  however, it is owned by exactly one of the processors
 *  on which it resides.
 */
class Entity {
public:

  /** \brief  The rank of this entity. */
  EntityRank entity_rank() const { return stk::mesh::entity_rank( m_key ); }

  /** \brief  The type (a.k.a. rank) of this entity. */
  EntityRank entity_type() const { return Entity::entity_rank(); }

  /** \brief  Identifier for this entity which is globally unique
   *          for a given entity type.
   */
  EntityId identifier() const { return stk::mesh::entity_id( m_key ); }

  /** \brief  The globally unique key ( entity type + identifier )
   *          of this entity.
   */
  const EntityKey & key() const { return m_key ; }

  /** \brief  The bucket which holds this mesh entity's field data */
  const Bucket & bucket() const { return *m_bucket ; }

  /** \brief  The ordinal for this entity within its bucket. */
  unsigned bucket_ordinal() const { return m_bucket_ord ; }

  /** \brief  The mesh bulk data synchronized_count when this entity's
   *          part membership was most recently modified.
   *
   *  If ( mesh.synchronized_state() == false &&
   *       mesh.synchronized_count() == entity.synchronized_count() )
   *  then entity was modified during this modification phase.
   */
  size_t synchronized_count() const { return m_sync_count ; }

  //------------------------------------
  /** \brief  All \ref stk::mesh::Relation "Entity relations"
   *          for which this entity is a member.
   */
  PairIterRelation relations() const { return PairIterRelation( m_relation ); }

  /** \brief  \ref stk::mesh::Relation "Entity relations" for which this
   *          entity is a member, the other entity is of a given type.
   */
  PairIterRelation relations( unsigned type ) const ;

  //------------------------------------
  /** \brief  Parallel processor rank of the processor which owns this entity */
  unsigned owner_rank() const { return m_owner_rank ; }

  /** \brief  Parallel processes which share this entity. */
  PairIterEntityComm sharing() const ;

  /** \brief  Complete communicaiton list for this entity */
  PairIterEntityComm comm() const { return PairIterEntityComm( m_comm ); }

  /** \brief  Subset communicaiton list for this entity */
  PairIterEntityComm comm( const Ghosting & ) const ;

  //------------------------------------

private:

  const EntityKey              m_key ;       ///< Globally unique key
  std::vector<Relation>        m_relation ;  ///< This entity's relationships
  std::vector<EntityCommInfo>  m_comm ;      ///< This entity's communications
  Bucket *                     m_bucket ;    ///< Bucket for the entity's field data
  unsigned              m_bucket_ord ;       ///< Ordinal within the bucket
  unsigned              m_owner_rank ;       ///< Owner processors' rank
  size_t                m_sync_count ;       ///< Last membership change

  ~Entity();
  explicit Entity( const EntityKey & arg_key );

  Entity(); ///< Default constructor not allowed
  Entity( const Entity & ); ///< Copy constructor not allowed
  Entity & operator = ( const Entity & ); ///< Assignment operator not allowed

  bool insert( const EntityCommInfo & );

#ifndef DOXYGEN_COMPILE
  friend class BulkData ;
#endif /* DOXYGEN_COMPILE */
};

/** \brief  Comparison operator for entities compares the entities' keys */
class EntityLess {
public:
  ~EntityLess() {}
  EntityLess() {}
  EntityLess( const EntityLess & ) {}
  EntityLess & operator = ( const EntityLess & ) { return *this ; }

  /** \brief  Comparison operator */
  bool operator()(const Entity& lhs, const Entity& rhs) const
  { return lhs.key() < rhs.key(); }

  bool operator()(const Entity& lhs, const EntityKey & rhs) const
  { return lhs.key() < rhs ; }

  /** \brief  Comparison operator */
  bool operator()(const Entity* lhs, const Entity* rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey() ;
    const EntityKey rhs_key = rhs ? rhs->key() : EntityKey() ;
    return lhs_key < rhs_key ;
  }

  bool operator()(const Entity* lhs, const Entity& rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey();
    return lhs_key < rhs.key() ;
  }

  bool operator()(const Entity* lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = lhs ? lhs->key() : EntityKey() ;
    return lhs_key < rhs ;
  }

  bool operator()( const EntityProc & lhs, const EntityProc & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    const EntityKey rhs_key = rhs.first ? rhs.first->key() : EntityKey() ;
    return lhs_key != rhs_key ? lhs_key < rhs_key : lhs.second < rhs.second ;
  }

  bool operator()( const EntityProc & lhs, const Entity & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    return lhs_key < rhs.key();
  }

  bool operator()( const EntityProc & lhs, const Entity * rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    const EntityKey rhs_key = rhs       ? rhs->key() : EntityKey() ;
    return lhs_key < rhs_key ;
  }

  bool operator()( const EntityProc & lhs, const EntityKey & rhs) const
  {
    const EntityKey lhs_key = lhs.first ? lhs.first->key() : EntityKey() ;
    return lhs_key < rhs ;
  }

}; //class EntityLess

class EntityEqual
{
public:
  bool operator()(const stk::mesh::Entity* lhs, const stk::mesh::Entity* rhs) const
  {
    const stk::mesh::EntityKey lhs_key = lhs ? lhs->key() : stk::mesh::EntityKey();
    const stk::mesh::EntityKey rhs_key = rhs ? rhs->key() : stk::mesh::EntityKey();
    return lhs_key == rhs_key;
  }

  bool operator()(const stk::mesh::Entity& lhs, const stk::mesh::Entity& rhs) const
  {
    const stk::mesh::EntityKey lhs_key = lhs.key();
    const stk::mesh::EntityKey rhs_key = rhs.key();
    return lhs_key == rhs_key;
  }
};

//----------------------------------------------------------------------

/** \} */

} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif /* stk_mesh_Entity_hpp */

