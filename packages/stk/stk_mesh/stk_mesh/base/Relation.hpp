#ifndef STK_MESH_RELATION_HPP
#define STK_MESH_RELATION_HPP

#include <stdint.h>                     // for uint32_t
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Types.hpp>      // for RelationType, EntityRank, etc
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowAssertMsg
#include <vector>                       // for vector
#include "boost/range/begin.hpp"        // for const_begin
#include "boost/range/end.hpp"          // for const_end
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */
//----------------------------------------------------------------------
/** \brief  A relation between two mesh entities with a
 *          relation <b> identifier </b> and <b> kind </b>.
 *
 *  Each entity owns a collection of relations to other entities.
 *  Each of these relations has a referenced entity, direction,
 *  relation identifier, and kind.
 *  - The direction specifies whether the relation is from the
 *    entity owning the relation to the referenced entity or
 *    conversely from the referenced entity to the owning entity.
 *  - The identifier provides a local numbering convention for relations,
 *    for example the local numbering convention for the nodes of an element.
 *  - Relations can be given a <em> kind </em> raw_relation_id to differentiate
 *    between topological relations (kind = 0) and other kinds of
 *    relationships such as constraints.
 *
 *  Relations are ordered by their
 *  - kind of relation,
 *  - type of the referenced entity,
 *  - if the relation is forward or converse,
 *  - local relation identifier, and then
 *  - referenced entity's identifier.
 */
class Relation {
public:
  typedef uint32_t raw_relation_id_type ;
  typedef uint32_t attribute_type;

  typedef stk::mesh::RelationType RelationType;

  /** \brief  Constructor */
  Relation();

  /** \brief  Construct a relation from a referenced entity and local identifier
   */
  Relation( const BulkData &mesh, Entity entity , RelationIdentifier identifier );
  // Defined in BulkData.hpp to break circular dependency.

  attribute_type   attribute() const { return m_attribute; }
  void set_attribute(attribute_type attr) const  { m_attribute = attr; }

  /** \brief  The encoded relation raw_relation_id */
  inline static raw_relation_id_type raw_relation_id( unsigned rank , unsigned id );

  /** \brief  The encoded relation raw_relation_id */
  raw_relation_id_type raw_relation_id() const { return m_raw_relation.value ; }

  /** \brief  The rank of the referenced entity */
  inline unsigned entity_rank() const ;

  /** \brief  The local relation identifier */
  inline RelationIdentifier relation_ordinal() const ;

  /** \brief  The referenced entity */
  inline Entity entity() const;

  /** \brief  Equality operator */
  inline bool operator == ( const Relation & r ) const;

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
  { return !(*this == r); }

private:

  enum {
    rank_digits = 8  ,
    id_digits   = 24 ,
    id_mask     = ~(0u) >> rank_digits
#ifdef SIERRA_MIGRATION
    ,
    fwmk_relation_type_digits = 8,
    fmwk_orientation_digits   = 24,
    fmwk_orientation_mask     = ~(0u) >> fwmk_relation_type_digits
#endif
  };

  union RawRelationType {
  public:
    raw_relation_id_type value ;

    struct {
      raw_relation_id_type identifier  : id_digits ;
      raw_relation_id_type entity_rank : rank_digits ;
    } normal_view ;

    struct {
      raw_relation_id_type entity_rank : rank_digits ;
      raw_relation_id_type identifier  : id_digits ;
    } reverse_view ;

    RawRelationType( raw_relation_id_type v ) : value(v) {}
    RawRelationType() : value(0) {}
    RawRelationType( const RawRelationType & rhs ) : value( rhs.value ) {}
    RawRelationType & operator = ( const RawRelationType & rhs )
      { value = rhs.value ; return *this ; }
  };

  RawRelationType        m_raw_relation ;
  mutable attribute_type m_attribute ;

  Entity                 m_target_entity;

 public:
  static raw_relation_id_type max_id() { return (1 << id_digits) - 1;}


// Issue: Framework supports relation types (parent, child, etc) that STK_Mesh
// does not support, so these relations will have to be managed by framework
// until:
// A) STK_Mesh is rewritten to support these extra relation types
// B) An extra data-structure is added to manage these relation types and
// this data structure works together with STK_mesh under the generic API
// to manage all relation-types.
// C) We transition to using a mesh that supports everything framework
// supports and this problem goes away.
//
// The problem is that, with framework managing some relations and STK_Mesh
// managing others, the type of the relation descriptor is different depending
// on what type of relation you're dealing with. This can be addressed with
// templates, but this makes the code very ugly. Instead...
//
// Solution: Have framework and STK_Mesh use the same type as its relation_descriptor.
// The code below is designed to make this class compatible with the fmwk
// Relation class.
#ifdef SIERRA_MIGRATION
 public:

  // Moved this to enum in struct RelationType.
  //   /**
  //   * Predefined identifiers for mesh object relationship types.
  //   */
  //  enum RelationType {
  //    USES	= 0 ,
  //    USED_BY	= 1 ,
  //    CHILD	= 2 ,
  //    PARENT	= 3 ,
  //    EMBEDDED	= 0x00ff , // 4
  //    CONTACT	= 0x00ff , // 5
  //    AUXILIARY   = 0x00ff ,
  //    INVALID     = 10
  //  };

  enum {
    POLARITY_MASK       = 0x80,
    POLARITY_POSITIVE   = 0x80,
    POLARITY_NEGATIVE   = 0x00,
    POLARITY_IDENTITY   = 0x80
  };

  static bool polarity(unsigned orient) {
    return (orient & POLARITY_MASK) == POLARITY_POSITIVE;
  }

  static unsigned permutation(unsigned orient) {
    return orient & ~POLARITY_MASK;
  }

  /**
   * Construct filled-out relation, fmwk-style
   */
  // Only needed by Framework and Framework-based apps.
  Relation(EntityRank, Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient = 0);

  // Only needed by Framework and Framework-based apps.
  inline void setMeshObj(Entity object, EntityRank object_rank);

  RelationType  getRelationType() const {
    return static_cast<RelationType::relation_type_t >(attribute() >> fmwk_orientation_digits);
  }

  void setRelationType(RelationType relation_type) {
    set_attribute( (relation_type << fmwk_orientation_digits) | getOrientation() );
  }

  RelationIdentifier getOrdinal() const {
    return relation_ordinal();
  }

  void setOrdinal(RelationIdentifier ordinal) {
    m_raw_relation = Relation::raw_relation_id( entity_rank(), ordinal );
  }

  attribute_type getOrientation() const {
    return attribute() & fmwk_orientation_mask;
  }

  void setOrientation(attribute_type orientation) {
    set_attribute( (getRelationType() << fmwk_orientation_digits) | orientation );
  }

  /**
   * Query polarity of the related mesh object.
   * A 'true' polarity indicates that the related mesh object is aligned.
   * For element-edge or face-edge a 'true' polarity indicates that the
   * nodes of the edge are compatibly ordered with the nodes of the
   * element or faces.
   * For element-face a 'true' polarity indicates that the ordering
   * of face-nodes is compatible with the ordering of element nodes,
   * i.e. the face's normal defined by a clockwise ordering is outward.
   */
  bool polarity() const
  { return compute_polarity(getOrientation()); }

  // TODO: This doesn't belong here. This is topology information.
  static bool compute_polarity(attribute_type orientation)
  { return (orientation & POLARITY_MASK) == POLARITY_POSITIVE; }

  unsigned permutation() const {
    return getOrientation() & ~POLARITY_MASK;
  }

private:
  bool has_fmwk_state() const { return getRelationType() != RelationType::INVALID; }
#endif // SIERRA_MIGRATION
};

//----------------------------------------------------------------------

inline
Relation::raw_relation_id_type
Relation::raw_relation_id( unsigned rank , unsigned id )
{
  ThrowAssertMsg( id <= id_mask,
                  "For args rank " << rank << ", id " << id << ": " <<
                  "id " << " > id_mask=" << id_mask );

  return ( raw_relation_id_type(rank) << id_digits ) | id ;
}

inline
unsigned Relation::entity_rank() const
{ return m_raw_relation.value >> id_digits; }

inline
RelationIdentifier Relation::relation_ordinal() const
{ return unsigned( m_raw_relation.value & id_mask ); }


//----------------------------------------------------------------------
/** \brief  Query which mesh entities have a relation
 *          to all of the input mesh entities.
 */
void get_entities_through_relations(
  const BulkData& mesh,
  const std::vector<Entity> & entities ,
        std::vector<Entity> & entities_related );

/** \brief  Query which mesh entities have a relation
 *          to all of the input mesh entities of the given
 *          mesh rank.
 */
void get_entities_through_relations(
  const BulkData& mesh,
  const std::vector<Entity> & entities ,
        EntityRank             entities_related_rank ,
        std::vector<Entity> & entities_related );

/** \brief  Induce entities' part membership based upon relationships
 *          between entities. Insert the result into 'induced_parts'.
 */
void induced_part_membership( const Part & part ,
                              EntityRank entity_rank_from ,
                              EntityRank entity_rank_to ,
                              OrdinalVector & induced_parts,
                              bool include_supersets=true);

/** \brief  Induce entities' part membership based upon relationships
 *          between entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership(const BulkData& mesh, const Entity entity_from ,
                              const OrdinalVector       & omit ,
                                    EntityRank            entity_rank_to ,
                                    OrdinalVector       & induced_parts,
                                    bool include_supersets=true);

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership(const BulkData& mesh, const Entity entity ,
                              const OrdinalVector & omit ,
                                    OrdinalVector & induced_parts,
                                    bool include_supersets=true);


//----------------------------------------------------------------------

#if 0
/** \brief  Decode and print the relation raw_relation_id */
std::ostream &
print_relation( std::ostream & , Relation::raw__attr_type );

/** \brief  Decode and print the relation raw_relation_id and referenced entity key */
std::ostream &
print_relation( std::ostream & , const MetaData & ,
                Relation::raw__attr_type , EntityKey );
#endif

/** \brief  Print the relation raw_relation_ids and referenced entity's key */
std::ostream & operator << ( std::ostream & , const Relation & );

/** \} */

inline
Relation::Relation() :
  m_raw_relation(),
  m_attribute(),
  m_target_entity()
{
#ifdef SIERRA_MIGRATION
  setRelationType(RelationType::INVALID);
#endif
}

inline
bool Relation::operator == ( const Relation & rhs ) const
{
  return m_raw_relation.value == rhs.m_raw_relation.value && m_target_entity == rhs.m_target_entity
#ifdef SIERRA_MIGRATION
    // compared fmwk state too
    && m_attribute == rhs.m_attribute
#endif
    ;
}

inline
bool same_specification(const Relation& lhs, const Relation& rhs)
{
#ifdef SIERRA_MIGRATION
  return  lhs.entity_rank()     == rhs.entity_rank() &&
          lhs.getRelationType() == rhs.getRelationType() &&
          lhs.getOrdinal()      == rhs.getOrdinal();
#else
  return  lhs.entity_rank()     == rhs.entity_rank();
#endif // SIERRA_MIGRATION
}

#ifdef SIERRA_MIGRATION

inline
RelationType
back_relation_type(const RelationType relType)
{
  /* %TRACE[NONE]% */  /* %TRACE% */
  switch(relType) {
  case RelationType::USES:
    return RelationType::USED_BY;
  case RelationType::USED_BY:
    return RelationType::USES;
  case RelationType::CHILD:
    return RelationType::PARENT;
  case RelationType::PARENT:
    return RelationType::CHILD;
  default:
    return relType;
  }
}

// Made publicly available so that MeshObj.C can use it
template <class Iterator>
bool
verify_relation_ordering(Iterator begin, Iterator end)
{
  for (Iterator itr = begin; itr != end; ) {
    Iterator prev = itr;
    ++itr;
    if (itr != end) {

      if (itr->entity_rank() < prev->entity_rank()) {
        return false ;
      }

      if (itr->entity_rank() == prev->entity_rank()) {

        if (itr->getRelationType() < prev->getRelationType()) {
          return false ;
        }

        if (itr->getRelationType() == prev->getRelationType()) {

          if (itr->getOrdinal() < prev->getOrdinal()) {
            return false ;
          }
        }
      }
    }
  }
  return true ;
}

template <class Range>
bool
verify_relation_ordering(Range range)
{
  return verify_relation_ordering(boost::const_begin(range), boost::const_end(range));
}

namespace impl {

inline
bool internal_is_handled_generically(const RelationType relation_type)
{ return relation_type == RelationType::USES || relation_type == RelationType::USED_BY; }

}

#endif // SIERRA_MIGRATION

inline
Entity Relation::entity() const
{
  return m_target_entity;
}


#ifdef SIERRA_MIGRATION

inline
Relation::Relation(EntityRank rel_rank, Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned orient)
  :
      m_raw_relation( Relation::raw_relation_id(rel_rank, ordinal )),
      m_attribute( (relation_type << fmwk_orientation_digits) | orient ),
      m_target_entity(obj)
{
  ThrowAssertMsg( orient <= fmwk_orientation_mask,
      "orientation " << orient << " exceeds maximum allowed value");
}

inline
void Relation::setMeshObj(Entity object, EntityRank object_rank )
{
  m_raw_relation = Relation::raw_relation_id( object_rank, relation_ordinal() );
  m_target_entity = object;
}

#endif

} // namespace mesh
} // namespace stk
#endif
