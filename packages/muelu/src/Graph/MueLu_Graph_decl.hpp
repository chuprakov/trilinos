#ifndef MUELU_GRAPH_DECL_HPP
#define MUELU_GRAPH_DECL_HPP

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_CrsGraph.hpp>     // inline functions requires class declaration
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"
#include "MueLu_Graph_fwd.hpp"

/******************************************************************************
   MueLu representation of a graph.
******************************************************************************/

namespace MueLu {

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class Graph 
    : public BaseClass {
#undef MUELU_GRAPH_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    Graph(const RCP<const CrsGraph> & graph, const std::string & objectLabel="") : graph_(graph) { 
      //setObjectLabel(objectLabel); 
      globalamalblockid2myrowid_ = Teuchos::null;
      globalamalblockid2globalrowid_ = Teuchos::null;
    }

    virtual ~Graph() {}
    
    size_t GetNodeNumVertices() const { return graph_->getNodeNumRows(); }
    size_t GetNodeNumEdges()    const { return graph_->getNodeNumEntries(); }
    
    Xpetra::global_size_t GetGlobalNumEdges() const { return graph_->getGlobalNumEntries(); }

    const RCP<const Teuchos::Comm<int> > GetComm() const { return graph_->getComm(); }
    const RCP<const Map> GetDomainMap() const { return graph_->getDomainMap(); }

    //! returns overlapping import map (nodes)
    const RCP<const Map> GetImportMap() const { return graph_->getColMap();    }

    //! Return the list of vertices adjacent to the vertex 'v'
    Teuchos::ArrayView<const LocalOrdinal> getNeighborVertices(LocalOrdinal v) const;

    //! store amalgamation information in MueLu::Graph object
    //! Both maps use the global block id of the amalgamated matrix as key and store the corresponding
    //! local row DOF ids and the global row DOF ids.
    //! The map globalamalblockid2globalrowid is used by MueLu::Graph::GetImportDofMap to generate
    //! the overlapping DofMap that is needed by the tentative prolongation operator (TentativePFactory).
    //! The map globalamalblockid2myrowid can be accessed from outside by GetAmalgamationParams and is needed
    //! by the aggregation algorithm
    void SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid,RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid) const;

    //! returns amalgamation information globalamalblockid2myrowid
    //! only valid if SetAmalgamationParams has been used before (i.e. CoalesceDropFactory::Amalgamate is called).
    RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > GetMyAmalgamationParams() const;

    //! returns amalgamation information globalamalblockid2globalrowid
    //! only valid if SetAmalgamationParams has been used before (i.e. CoalesceDropFactory::Amalgamate is called).
    RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > GetGlobalAmalgamationParams() const;

#ifdef MUELU_UNUSED
    size_t GetNodeNumGhost() const;
#endif

    /// Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

  private:

    RCP<const CrsGraph> graph_;

    //! @name amalgamation information variables
    //@{

    /// map: global block id of amalagamated matrix -> vector of local row ids of unamalgamated matrix (only for global block ids of current proc)
    mutable RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid_;   //< used by aggregation factory
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid_; //< used for building overlapping ImportDofMap

    //@}

  };

} // namespace MueLu

#define MUELU_GRAPH_SHORT
#endif // MUELU_GRAPH_DECL_HPP
