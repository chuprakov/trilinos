/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/Selector.hpp>
#include <iostream>                     // for operator<<
#include <stk_mesh/base/Bucket.hpp>     // for has_superset
#include <stk_mesh/base/Types.hpp>      // for PartVector
#include <string>                       // for operator<<
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/FieldRestriction.hpp"
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowRequireMsg


namespace stk {
namespace mesh {

namespace {

using impl::SelectorNode;

const char* to_str(SelectorNodeType::node_type type)
{
  switch(type) {
  case SelectorNodeType::UNION:
    return " | ";
  case SelectorNodeType::INTERSECTION:
    return " & ";
  default:
    return " - ";
  };
}

std::ostream& print_expr_impl(std::ostream & out, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
  case SelectorNodeType::INTERSECTION:
  case SelectorNodeType::DIFFERENCE:
    out << "(";
    print_expr_impl(out, root->lhs());
    out << to_str(root->m_type);
    print_expr_impl(out, root->rhs());
    out << ")";
    break;
  case SelectorNodeType::COMPLEMENT:
    out << "!(";
    print_expr_impl(out, root->unary());
    out << ")";
    break;
  case SelectorNodeType::PART:
    if (root->part() != NULL) {
      out << root->part()->name();
    }
    else {
      out << "NOTHING";
    }
    break;
  };
  return out;
}

bool select_bucket_impl(Bucket const& bucket, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_bucket_impl(bucket, root->lhs()) || select_bucket_impl(bucket, root->rhs());
  case SelectorNodeType::INTERSECTION:
    return select_bucket_impl(bucket, root->lhs()) && select_bucket_impl(bucket, root->rhs());
  case SelectorNodeType::DIFFERENCE:
    return select_bucket_impl(bucket, root->lhs()) && !select_bucket_impl(bucket, root->rhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_bucket_impl(bucket, root->unary());
  case SelectorNodeType::PART:
    return (root->part() != NULL)? has_superset(bucket, *root->part()) : false;
  default:
    return false;
  };
}

bool select_part_impl(Part const& part, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_part_impl(part, root->lhs()) || select_part_impl(part, root->rhs());
  case SelectorNodeType::INTERSECTION:
    return select_part_impl(part, root->lhs()) && select_part_impl(part, root->rhs());
  case SelectorNodeType::DIFFERENCE:
    return select_part_impl(part, root->lhs()) && !select_part_impl(part, root->rhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_part_impl(part, root->unary());
  case SelectorNodeType::PART:
    return (root->part() != NULL) ? root->part()->contains(part) : false;
  default:
    return false;
  };
}

bool is_all_union_impl(SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return is_all_union_impl(root->lhs()) && is_all_union_impl(root->rhs());
  case SelectorNodeType::INTERSECTION:
  case SelectorNodeType::DIFFERENCE:
  case SelectorNodeType::COMPLEMENT:
    return false;
  case SelectorNodeType::PART:
    return root->part() != NULL;
  default:
    return false;
  };
}

void gather_parts_impl(PartVector& parts, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    gather_parts_impl(parts, root->lhs());
    gather_parts_impl(parts, root->rhs());
    break;
  case SelectorNodeType::INTERSECTION:
    // HACK: Only first part (picks up Context part)
    gather_parts_impl(parts, root->lhs());
    break;
  case SelectorNodeType::DIFFERENCE:
    ThrowRequireMsg(false, "Cannot get_parts from a selector with differences");
    break;
  case SelectorNodeType::COMPLEMENT:
    ThrowRequireMsg(false, "Cannot get_parts from a selector with differences");
    break;
  case SelectorNodeType::PART:
    if (root->part() != NULL) parts.push_back(const_cast<Part*>(root->part()));
  };
}

bool select_part_vector_impl(PartVector const& parts, SelectorNode const* root)
{
  switch(root->m_type) {
  case SelectorNodeType::UNION:
    return select_part_vector_impl(parts, root->lhs()) || select_part_vector_impl(parts, root->rhs());
  case SelectorNodeType::INTERSECTION:
    return select_part_vector_impl(parts, root->lhs()) && select_part_vector_impl(parts, root->rhs());
  case SelectorNodeType::DIFFERENCE:
    return select_part_vector_impl(parts, root->lhs()) && !select_part_vector_impl(parts, root->rhs());
  case SelectorNodeType::COMPLEMENT:
    return !select_part_vector_impl(parts, root->unary());
  case SelectorNodeType::PART:
    if (root->part() == NULL) {
      return false;
    }
    else {
      for (size_t i = 0, ie = parts.size(); i < ie; ++i) {
        if (parts[i] == root->part()) {
          return true;
        }
      }
      return false;
    }
  default:
    return false;
  };
}

} // namespace

std::ostream & operator<<( std::ostream & out, const Selector & selector)
{
  return print_expr_impl(out, &selector.m_expr[0]);
}


bool Selector::operator()( const Part & part ) const
{
  return select_part_impl(part, &m_expr[0]);
}

bool Selector::operator()( const Part * part ) const
{
  return select_part_impl(*part, &m_expr[0]);
}

bool Selector::operator()( const Bucket & bucket ) const
{
  return select_bucket_impl(bucket, &m_expr[0]);
}

bool Selector::operator()( const Bucket * bucket ) const
{
  return select_bucket_impl(*bucket, &m_expr[0]);
}

bool Selector::operator()(const PartVector& parts) const
{
  return select_part_vector_impl(parts, &m_expr[0]);
}

bool Selector::operator<(const Selector& rhs) const
{
  // kinda arbitrary, but should work as long as all we need is a consistent ordering
  if (m_expr.size() != rhs.m_expr.size()) {
    return m_expr.size() < rhs.m_expr.size();
  }

  for (size_t i = 0, ie = m_expr.size(); i < ie; ++i) {
    if (m_expr[i].m_type != rhs.m_expr[i].m_type) {
      return m_expr[i].m_type < rhs.m_expr[i].m_type;
    }
    if (m_expr[i].m_type == SelectorNodeType::PART &&
        m_expr[i].part() != rhs.m_expr[i].part()) {
      Part const* lhs_part = m_expr[i].part();
      Part const* rhs_part = rhs.m_expr[i].part();

      if (lhs_part != NULL && rhs_part != NULL) {
        return lhs_part->mesh_meta_data_ordinal() < rhs_part->mesh_meta_data_ordinal();
      }
      else if (lhs_part == NULL && rhs_part != NULL) {
        return false;
      }
      else if (lhs_part != NULL && rhs_part == NULL) {
        return true;
      }
    }
  }

  return false;
}

void Selector::get_parts(PartVector& parts) const
{
  gather_parts_impl(parts, &m_expr[0]);
}

BulkData* Selector::find_mesh() const
{
    BulkData* mesh = NULL;
    for(size_t i=0; i<m_expr.size(); ++i) {
        if (m_expr[i].node_type() == SelectorNodeType::PART && m_expr[i].part() != NULL) {
            mesh = &m_expr[i].part()->mesh_bulk_data();
        }
    }
    return mesh;
}

BucketVector const& Selector::get_buckets(EntityRank entity_rank) const
{
    static BucketVector emptyBucketVector;
    if (m_expr.empty()) {
        return emptyBucketVector;
    }

    BulkData* mesh = find_mesh();
    ThrowRequireMsg(mesh != NULL,
        "ERROR, Selector::get_buckets not available if selector expression does not involve any mesh Parts.");

    return mesh->get_buckets(entity_rank, *this);
}

bool Selector::is_empty(EntityRank entity_rank) const
{
    if (m_expr.empty()) {
        return true;
    }

    BulkData * mesh = this->find_mesh();
    ThrowRequireMsg(mesh != NULL,
                    "ERROR, Selector::empty not available if selector expression does not involve any mesh Parts.");
    if (mesh->synchronized_state() == BulkData::MODIFIABLE) {
      BucketVector const& buckets = this->get_buckets(entity_rank);
      for(size_t i=0; i<buckets.size(); ++i) {
          if (buckets[i]->size() >0) {
              return false;
          }
      }
      return true;
    }
    return get_buckets(entity_rank).empty();
}


bool Selector::is_all_unions() const
{
  return is_all_union_impl(&m_expr[0]);
}

Selector selectUnion( const PartVector& union_part_vector )
{
  Selector selector;
  if (union_part_vector.size() > 0) {
    selector = *union_part_vector[0];
    for (unsigned i = 1 ; i < union_part_vector.size() ; ++i) {
      selector |= *union_part_vector[i];
    }
  }
  return selector;
}

Selector selectIntersection( const PartVector& intersection_part_vector )
{
  Selector selector;
  if (intersection_part_vector.size() > 0) {
    selector = *intersection_part_vector[0];
    for (unsigned i = 1 ; i < intersection_part_vector.size() ; ++i) {
      selector &= *intersection_part_vector[i];
    }
  }
  return selector;
}

Selector selectField( const FieldBase& field )
{
  Selector selector;
  const FieldRestrictionVector& sel_rvec = field.restrictions();
  for(size_t i=0; i<sel_rvec.size(); ++i) {
    selector |= sel_rvec[i].selector();
  }

  return selector;
}

bool is_subset(Selector const& lhs, Selector const& rhs)
{
  // If either selector has complements or intersections, it becomes
  // much harder to determine if one is a subset of the other
  if (lhs.is_all_unions() && rhs.is_all_unions()) {
    PartVector lhs_parts, rhs_parts;
    lhs.get_parts(lhs_parts);
    rhs.get_parts(rhs_parts);
    for (size_t l = 0, le = lhs_parts.size(); l < le; ++l) {
      Part const& lhs_part = *lhs_parts[l];
      bool found = false;
      for (size_t r = 0, re = rhs_parts.size(); !found && r < re; ++r) {
        Part const& rhs_part = *rhs_parts[r];
        found = rhs_part.contains(lhs_part);
      }
      if (!found) {
        return false;
      }
    }
    return true;
  }
  else {
    return false;
  }
}

} // namespace mesh
} // namespace stk
