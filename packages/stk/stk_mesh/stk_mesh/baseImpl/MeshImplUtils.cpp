/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

//----------------------------------------------------------------------

#include <stk_mesh/base/BulkData.hpp>
#include <stk_topology/topology.hpp>

#include <vector>

//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

void find_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems)
{
  elems.clear();
  std::vector<Entity> tmp;
  for(unsigned i=0; i<numNodes; ++i) {
    const Entity* elements = mesh.begin_elements(nodes[i]);
    unsigned numElements = mesh.num_elements(nodes[i]);
    tmp.assign(elements, elements+numElements);
    std::sort(tmp.begin(), tmp.end());
    if (i==0) {
      elems.assign(tmp.begin(), tmp.end());
    }
    else {
       std::vector<Entity> intersect;
       std::back_insert_iterator<std::vector<Entity> > intersect_itr(intersect);
       std::set_intersection(elems.begin(), elems.end(),
                             tmp.begin(), tmp.end(),
                             intersect_itr);
       elems.swap(intersect);
    }
  }
}

void find_locally_owned_elements_these_nodes_have_in_common(BulkData& mesh, unsigned numNodes, const Entity* nodes, std::vector<Entity>& elems)
{
  find_elements_these_nodes_have_in_common(mesh, numNodes, nodes, elems);

  for(int i=elems.size()-1; i>=0; --i) {
    if (!mesh.bucket(elems[i]).owned()) {
      elems.erase(elems.begin()+i);
    }
  }
}

bool find_element_edge_ordinal_and_equivalent_nodes(BulkData& mesh, Entity element, unsigned numEdgeNodes, const Entity* edgeNodes, unsigned& elemEdgeOrdinal, Entity* elemEdgeNodes)
{
  stk::topology elemTopology = mesh.bucket(element).topology();
  stk::topology edgeTopology = elemTopology.edge_topology();
  const Entity* elemNodes = mesh.begin_nodes(element);
  ThrowAssertMsg(mesh.num_nodes(element) == elemTopology.num_nodes(), "findElementEdgeOrdinalAndNodes ERROR, element (id="<<mesh.identifier(element)<<") has wrong number of connected nodes ("<<mesh.num_nodes(element)<<"), expected elemTopology.num_nodes()="<<elemTopology.num_nodes());

  unsigned numEdgesPerElem = elemTopology.num_edges();
  for(elemEdgeOrdinal=0; elemEdgeOrdinal<numEdgesPerElem; ++elemEdgeOrdinal) {
    elemTopology.edge_nodes(elemNodes, elemEdgeOrdinal, elemEdgeNodes);
    if (edgeTopology.equivalent(edgeNodes, elemEdgeNodes).first) {
      //found element edge equivalent to edgeNodes.
      //output arguments elemEdgeOrdinal and elemEdgeNodes are set, let's get out of here.
      return true;
    }
  }

  return false;//didn't find element edge equivalent to input edgeNodes
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

