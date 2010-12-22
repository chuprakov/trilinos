#ifndef MUELU_GRAPH_HPP
#define MUELU_GRAPH_HPP

#include <Teuchos_ArrayView.hpp>
#include <Cthulhu_CrsGraph.hpp>

#include "MueLu_Exceptions.hpp"

/******************************************************************************
   MueLu representation of a graph.
******************************************************************************/

namespace MueLu {
  
  template <class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps >
  class Graph 
    : public Teuchos::Describable {
    
#include "MueLu_UseShortNames.hpp"

  public:

    Graph(const RCP<const CrsGraph> & graph, const std::string & objectLabel="") : graph_(graph) { setObjectLabel(objectLabel); }
    ~Graph() {}
    
    inline int GetNodeNumVertices() const { return graph_->getNodeNumRows(); }
    inline int GetNodeNumEdges() const    { return graph_->getNodeNumEntries();  }
    
    inline int GetGlobalNumEdges() const { return graph_->getGlobalNumEntries(); }

    inline const RCP<const Teuchos::Comm<int> > GetComm() const { return graph_->getComm(); }
    inline const RCP<const Map> GetDomainMap() const { return graph_->getDomainMap(); }
    inline const RCP<const Map> GetImportMap() const { return graph_->getColMap(); }

    //! Return the list of vertices adjacent to the vertex 'v'
    inline Teuchos::ArrayView<const int> getNeighborVertices(int v) const { 
      Teuchos::ArrayView<const int> neighborVertices;
      graph_->getLocalRowView(v, neighborVertices); 
      return neighborVertices;
    }

    int GetNodeNumGhost() const { 
      /*
        Ray's comments about nGhost:
        Graph->NGhost == graph_->RowMatrixColMap()->NumMyElements() - graph_->OperatorDomainMap()->NumMyElements()
        is basically right. But we've had some issues about how epetra handles empty columns.
        Probably worth discussing this with Jonathan and Chris to see if this is ALWAYS right. 
      */
      int nGhost;//TODO = graph_->getColMap()->getNodeNumElements() - graph_->getDomainMap()->getNodeNumElements();
      if (nGhost < 0) nGhost = 0;
      
      return nGhost;
    }

  private:

    RCP<const CrsGraph> graph_;

  };

}

#define MUELU_GRAPH_SHORT
#endif
