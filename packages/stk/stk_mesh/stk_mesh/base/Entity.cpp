/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <stddef.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <boost/mpl/assert.hpp>

#ifdef SIERRA_MIGRATION
namespace {
static const std::vector<stk::mesh::Relation> dummy_vector;
}

namespace sierra {
namespace Fmwk {

const unsigned int INVALID_LOCAL_ID = std::numeric_limits<unsigned int>::max();
const stk::mesh::RelationIterator INVALID_RELATION_ITR = dummy_vector.end(); // Some STL implementation use POD for iterators

unsigned get_derived_type(const stk::mesh::Entity );

}
}
#endif

namespace stk {
namespace mesh {

// TODO - Activate once we move to intel-12.1
//BOOST_MPL_ASSERT(( boost::is_pod<Entity> ));

//----------------------------------------------------------------------

std::string print_entity_key(const Entity entity)
{
  return print_entity_key(MetaData::get(entity),
                          entity.key());
}

//
//----------------------------------------------------------------------

#ifdef SIERRA_MIGRATION

std::string Entity::TypeToString (Entity::ObjectTypeEnum type)
{
  if(type == NODE      ) return "NODE";
  if(type == EDGE      ) return "EDGE";
  if(type == FACE      ) return "FACE";
  if(type == ELEMENT   ) return "ELEMENT";
  if(type == CONSTRAINT) return "CONSTRANT";
  if(type == BASE_CLASS) return "BASE_CLASS";
  return "UNKNOWN";
}

// ---------------------------------------------------------------------

void Entity::internal_swap_in_real_entity(const int globalId)
{
  ThrowRequire(globalId > 0);
  m_entityImpl->m_fmwk_attrs.global_id  = globalId;

  BulkData::get(*this).change_entity_id(globalId, *this);

  internal_verify_initialization_invariant();

  // Issue: Fmwk-managed relations (also called auxiliary relations, are not
  // being resorted here, so we have to use a different < operator
  // when looking for relations

#ifndef NDEBUG
  internal_verify_meshobj_invariant();
#endif
}

// ---------------------------------------------------------------------

namespace {

struct IgnoreIdOrder
{
  bool operator()(const Relation& lhs, const Relation& rhs)
  {
    bool result = false;

    if (lhs.entity_rank() != rhs.entity_rank()) {
      result = lhs.entity_rank() < rhs.entity_rank();
    }
    else if (lhs.getRelationType() != rhs.getRelationType()) {
      result = lhs.getRelationType() < rhs.getRelationType();
    }
    else {
      result = lhs.relation_ordinal() < rhs.relation_ordinal();
    }
    return result;
  }
};

}

RelationIterator Entity::find_relation(const Relation& relation) const
{
  // Extremely hacky: It would be better to set up the < operator for relations so that lower_bound
  // can return the desired iterator, but any sane definition would probably force a change in
  // relation ordering and that's more than I'm willing to take on now.
  //
  // The current semantics for relation-searching is as follows:
  // Ordered based on derived_type, relation_type, and ordinal in descending precedence
  //   If multiple relations have the same derived_type, relation_type, and ordinal, a linear
  //   scan takes place looking for a matching meshobj. If no such meshobj was found, then
  //   we are left with an iterator pointing to the first relation with a different derived_type,
  //   relation_type, or ordinal. To sum up, the result of the search can either be equivalent to
  //   lower_bound OR upper_bound depending upon the state of the relations... YUCK!

  const Relation::RelationType relation_type = relation.getRelationType();

  RelationIterator rel = std::lower_bound(internal_begin_relation(relation_type),
                                          internal_end_relation(relation_type),
                                          relation,
                                          IgnoreIdOrder());

  // Should only loop if we are looking at back-relations, otherwise, relations with
  // matching specifications are not legal.
  while (rel != internal_end_relation(relation_type) &&
         same_specification(*rel, relation) &&
         rel->entity()      != relation.entity())
    ++rel;

  return rel;
}

// ---------------------------------------------------------------------

bool Entity::update_relation(
  const RelationIterator        ir ,
  const bool                    back_rel_flag) const
{
  const Relation::RelationType relType = ir->getRelationType();
  ThrowAssert(verify_relation_ordering(internal_begin_relation(relType), internal_end_relation(relType)));
  ThrowAssertMsg(!internal_is_handled_generically(relType),
                 "update_relation should not be called for STK-managed relations");

  Entity meshObj = ir->entity();

  const Relation::RelationType backRelType = back_relation_type(relType);

  ThrowAssert(verify_relation_ordering(meshObj.internal_begin_relation(backRelType), meshObj.internal_end_relation(backRelType)));

  // Create the corresponding back relation to ir
  Relation backRel_obj(*this, backRelType, ir->getOrdinal(), ir->getOrientation());
  RelationIterator backRel_itr = meshObj.find_relation(backRel_obj);

  const bool exists = backRel_itr != meshObj.internal_end_relation(backRelType) && *backRel_itr == backRel_obj;

  if (exists && !back_rel_flag) {
    // Remove back relation and increment the counter

    meshObj.erase_and_clear_if_empty(backRel_itr);

    //ThrowAssert(sierra::Fmwk::get_derived_type(meshObj) != Entity::ELEMENT);

    meshObj.inc_connection();
  }
  else if (!exists && back_rel_flag) {
    // Insert back relation

    const unsigned k = backRel_itr - meshObj.internal_begin_relation(backRelType) ;

    meshObj.reserve_relation(meshObj.aux_relations().size() + 1);

    meshObj.aux_relations().insert(meshObj.aux_relations().begin() + k, backRel_obj);

    //ThrowAssert(sierra::Fmwk::get_derived_type(meshObj) != Entity::ELEMENT);

    meshObj.dec_connection();
  }

  ThrowAssert(verify_relation_ordering(meshObj.internal_begin_relation(relType), meshObj.internal_end_relation(relType)));

  return true;
}

// ---------------------------------------------------------------------

void Entity::internal_verify_meshobj_invariant() const
{
  PairIterRelation stk_relations = relations();
  for ( ; !stk_relations.empty(); ++stk_relations ) {
    ThrowRequireMsg(stk_relations->entity().is_valid(), "Problem with: " << *stk_relations);
  }

  RelationVector& aux_relations = m_entityImpl->m_fmwk_attrs.aux_relations;
  for (RelationVector::const_iterator itr = aux_relations.begin(), end = aux_relations.end(); itr != end; ++itr) {
    ThrowRequireMsg(itr->entity().is_valid(), "Problem with: " << *itr);
  }
}

// ---------------------------------------------------------------------

void Entity::set_relation_orientation(RelationIterator rel, unsigned orientation)
{
  const Relation::RelationType backRelType = back_relation_type(rel->getRelationType());

  Entity meshObj = rel->entity();
  Relation backRel_obj(*this, backRelType, rel->getOrdinal(), rel->getOrientation());
  RelationIterator backRel_itr = meshObj.find_relation(backRel_obj);

  const bool exists = backRel_itr != meshObj.internal_end_relation(backRelType) && *backRel_itr == backRel_obj;
  ThrowRequire(exists);

  // Allow clients to make changes to orientation
  // Orientations do not affect Relation ordering, so this is safe.
  const_cast<Relation*>(&*rel)->setOrientation(orientation);
  const_cast<Relation*>(&*backRel_itr)->setOrientation(orientation);
}

#endif

BOOST_STATIC_ASSERT(( (int)MetaData::NODE_RANK == (int)Entity::NODE ));
BOOST_STATIC_ASSERT(( (int)MetaData::EDGE_RANK == (int)Entity::EDGE ));
BOOST_STATIC_ASSERT(( (int)MetaData::FACE_RANK == (int)Entity::FACE ));
BOOST_STATIC_ASSERT(( (int)MetaData::ELEMENT_RANK == (int)Entity::ELEMENT ));

} // namespace mesh
} // namespace stk
