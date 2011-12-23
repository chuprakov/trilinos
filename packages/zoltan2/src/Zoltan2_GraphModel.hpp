// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_GraphModel.hpp

    \brief The interface and implementations of a graph model.
*/


#ifndef _ZOLTAN2_GRAPHMODEL_HPP_
#define _ZOLTAN2_GRAPHMODEL_HPP_

#include <Zoltan2_Model.hpp>
#include <Zoltan2_MatrixInput.hpp>

#include <vector>
#include <Teuchos_Hashtable.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphModel
    \brief GraphModel defines the interface required for graph models.  

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    methods should be local methods.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  
  GraphModel(){
    throw std::logic_error("in non-specialized GraphModel");
  }

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return 0; }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const { return 0; }

  /*! Returns the number of global edges on this process.
   *  Includes remote edges.
   */
  size_t getLocalNumGlobalEdges() const { return 0; }

  /*! Returns the number of local edges on this process.
   *  Does not include remote edges.
   */
  size_t getLocalNumLocalEdges() const { return 0; }

  /*! Returns the global number edges.
   */
  global_size_t getGlobalNumEdges() const { return 0; }

  /*! Returns the dimension (0 or greater) of vertex weights.
   */
  int getVertexWeightDim() const { return 0; }

  /*! Returns the dimension (0 or greater) of edge weights.
   */
  int getEdgeWeightDim() const { return 0; }

  /*! Returns the dimension (0 to 3) of vertex coordinates.
   */
  int getCoordinateDim() const { return 0; }

  /*! Sets pointers to this process' vertex Ids and their weights.
      \param Ids will on return point to the list of the global Ids for
        each vertex on this process.
      \param xyz will on return point to a list coordinates for
         each vertex in the Ids list.  Coordinates are listed by
         vertex by component.
      \param wgts will on return point to a list of the weight or weights
         associated with each vertex in the Ids list.  Weights are listed by
         vertex by weight component.
       \return The number of ids in the Ids list.
      KDDKDD This function should not return coordinates, as coordinates are
      KDDKDD not needed in most graph algorithms.  There should be a separate
      KDDKDD function for coordinates.  Or there should be an option to 
      KDDKDD specify whether or not coordinates are needed.
   */
  // GNOs are returned.
  // LNOs are implied by the order that the vertices have been returned.

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const {
      return 0; }

  /*! getEdgeList:
      Sets pointers to this process' edge (neighbor) global Ids, including
      off-process edges.
      \param edgeIds This is the list of global neighbor Ids corresponding
        to the vertices listed in getVertexList.
      \param procIds lists the process owning each neighbor in the edgeIds
         list.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex.
      \param wgts will on return point to a list of the weight or weights
         associated with each edge in the edgeIds list.  Weights are listed by
         edge by weight component.
       \return The number of ids in the edgeIds list.
   */
  // Implied Vertex LNOs from getVertexList are used as indices to offsets 
  // array.
  // Vertex GNOs are returned as neighbors in edgeIds.

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  /*! getLocalEdgeList:
      Sets pointers to this process' local-only edge (neighbor) LNOs, using
      the same implied vertex LNOs returned in getVertexList.
      \param edgeIds This is the list of local neighbor Ids corresponding
        to the vertices listed in getVertexList.
        Off-process edges are not returned.
      \param offsets offsets[i] is the offset into edgeIds to the start
        of neighbors for ith vertex.
      \param wgts will on return point to a list of the weight or weights
         associated with each edge in the edgeIds list.  Weights are listed by
         edge by weight component.
       \return The number of ids in the edgeIds list.
   */
  // Implied Vertex LNOs from getVertexList are used as indices to offsets 
  // array.
  // Vertex LNOs are returned as neighbors in edgeIds.

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumVertices();
  }

  global_size_t getGlobalNumObjects() const
  {
    return getGlobalNumVertices();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const { return; }

  int getNumWeights() const { return 0; }
};

////////////////////////////////////////////////////////////////
// Graph model derived from MatrixInput.
//
//   TODO: support a flag that says the vertices are columns or
//           non-zeros rather than rows
////////////////////////////////////////////////////////////////

/*! Zoltan2::GraphModel<MatrixInput>
    \brief A specialization of GraphModel for Zoltan2::MatrixInput.
*/

template <typename User>
class GraphModel<MatrixInput<User> > : public Model<MatrixInput<User> >
{
public:

  typedef typename MatrixInput<User>::scalar_t  scalar_t;
  typedef typename MatrixInput<User>::gno_t     gno_t;
  typedef typename MatrixInput<User>::lno_t     lno_t;
  typedef typename MatrixInput<User>::gid_t     gid_t;
  typedef typename MatrixInput<User>::node_t    node_t;
  typedef IdentifierMap<User>     idmap_t;

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   *  \param  inputAdapter  an encapsulation of the user data
   *  \param  env           environment (library configuration settings)
   *  \param  consecutiveIdsRequired  set to true if the algorithm or
   *           third party library requires consecutive global vertex Ids.
   *  \param removeSelfEdges set to true if the algorithm or the third party
   *           library cannot handle self edges
   */
  GraphModel(const MatrixInput<User> *ia,
    const RCP<const Environment> &env, bool consecutiveIdsRequired=false,
    bool removeSelfEdges=false) :
     input_(ia), env_(env), gids_(), gnos_(), edgeGnos_(), procIds_(), 
     offsets_(), gnosConst_(), edgeGnosConst_(), procIdsConst_(), 
     numLocalEdges_(0), numGlobalEdges_(0), numLocalVtx_(0), 
     gidsAreGnos_(false), nearEdgeLnos_(), nearEdgeOffsets_(), 
     numNearLocalEdges_(0)

  {
    // Get the matrix from the input adapter

    gid_t const *vtxIds=NULL, *nborIds=NULL;
    lno_t const  *offsets=NULL;
    try{
      numLocalVtx_ = input_->getRowListView(vtxIds, offsets, nborIds);
    }
    Z2_FORWARD_EXCEPTIONS;

    gids_ = arcp(vtxIds, 0, numLocalVtx_, false);

    numLocalEdges_ = offsets[numLocalVtx_];

    // If Matrix has diagonal entries, and self-edges are to be removed
    //    do that now.

    ArrayRCP<lno_t> tmpOffsets;
    ArrayRCP<gid_t> tmpEdges;
    lno_t nSelfEdges = 0;

    size_t numOffsets = numLocalVtx_ + 1;

    if (removeSelfEdges) {

      lno_t *offArray = new lno_t [numOffsets];
      Z2_LOCAL_MEMORY_ASSERTION(*env, numOffsets, offArray);
      gid_t *edArray = new gid_t [numLocalEdges_];
      Z2_LOCAL_MEMORY_ASSERTION(*env, numLocalEdges_, !numLocalEdges_||edArray);

      for (lno_t i=0; i < numLocalVtx_; i++){

        offArray[i] = offsets[i] - nSelfEdges;

        for (lno_t j = offsets[i]; j < offsets[i+1]; j++) {
          if (gids_[i] == nborIds[j]) { // self edge; remove it
            nSelfEdges++;
          }
          else {  // Not a self-edge; keep it.
            edArray[j-nSelfEdges] = nborIds[j];
          }
        }
      }
      numLocalEdges_ -= nSelfEdges;
      offArray[numLocalVtx_] = numLocalEdges_;

      if (nSelfEdges > 0){
        tmpOffsets = arcp(offArray, 0, numLocalVtx_+1);
        tmpEdges = arcp(edArray, 0, numLocalEdges_);
      }
      else{
        delete [] offArray;
        if (numLocalEdges_) delete [] edArray;
      }
    }

    if (nSelfEdges == 0){
      offsets_ = arcp(const_cast<lno_t *>(offsets), 0, numOffsets, false);
      edgeGids_ = arcp(const_cast<gno_t *>(nborIds), 0, numLocalEdges_, false);
    }
    else{
      offsets_ = tmpOffsets;
      edgeGids_ =  tmpEdges;
    }

    reduceAll<int, size_t>(*(env_->comm_), Teuchos::REDUCE_SUM, 1,
      &numLocalEdges_, &numGlobalEdges_);

    // Create an IdentifierMap, which will map the user's global IDs to
    //   Zoltan2 internal global numbers if neccesary.
    //   The map can also give us owners of our vertex neighbors.

    RCP<const idmap_t> idMap;

    try{
      idMap = rcp(new idmap_t(env, gids_, consecutiveIdsRequired));
    }
    Z2_FORWARD_EXCEPTIONS;

    gidsAreGnos_ = idMap->gnosAreGids();

    if (numLocalVtx_ && !gidsAreGnos_){
      gno_t *tmp = new gno_t [numLocalVtx_];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalVtx_, tmp)
      gnos_ = arcp(tmp, 0, numLocalVtx_);

      try{
        // Because gidTranslate can translate gids to gnos or
        // gnos to gids depending on a flag, neither the gids nor
        // the gnos are declared to be const.
        ArrayRCP<gid_t> tmpGids = arcp_const_cast<gid_t>(gids_);

        idMap->gidTranslate(tmpGids(0,numLocalVtx_), gnos_(0,numLocalVtx_), 
          TRANSLATE_APP_TO_LIB);
      }
      Z2_FORWARD_EXCEPTIONS;

      if (numLocalEdges_){
        tmp = new gno_t [numLocalEdges_];
        Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalEdges_, tmp)
        edgeGnos_ = arcp(tmp, 0, numLocalEdges_);
      }
    }

    if (numLocalEdges_){
      int *p = new int [numLocalEdges_];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, numLocalEdges_, p)
      procIds_ = arcp(p, 0, numLocalEdges_);
    }

    if (!gidsAreGnos_){
      try{
        // All processes must make this call.
        idMap->gidGlobalTranslate(edgeGids_(0,numLocalEdges_), 
          edgeGnos_(0,numLocalEdges_), procIds_(0,numLocalEdges_));
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    this->setIdentifierMap(idMap);   // Zoltan2::Model method

    gnosConst_ = arcp_const_cast<const gno_t>(gnos_);
    edgeGnosConst_ = arcp_const_cast<const gno_t>(edgeGnos_);
    procIdsConst_ = arcp_const_cast<const int>(procIds_);
  }

  //!  Destructor
  ~GraphModel() { }

  // // // // // // // // // // // // // // // // // // // // // /
  // The GraphModel interface.
  // // // // // // // // // // // // // // // // // // // // // /

  size_t getLocalNumVertices() const
  {
    return input_->getLocalNumRows();
  }

  global_size_t getGlobalNumVertices() const
  {
    return input_->getGlobalNumRows();
  }

  size_t getLocalNumEdges() const
  {
    return numLocalEdges_;
  }

  global_size_t getGlobalNumEdges() const
  {
    return numGlobalEdges_;
  }

  int getVertexWeightDim() const
  {
    return 0;   // TODO
  }

  int getEdgeWeightDim() const
  {
    return 0;   // TODO
  }

  int getCoordinateDim() const
  {
    return 0;   // TODO
  }

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const
  {
    size_t n = getLocalNumVertices();
    Ids = ArrayView<const gno_t>(Teuchos::null);

    if (n){
      if (gidsAreGnos_)
        Ids = gids_(0, n);
      else
        Ids = gnosConst_(0, n);
    }

    return n;
    // KDDKDD  Is it dangerous that xyz is ignored here?  Perhaps coordinates
    // KDDKDD  should be separate.
    // LRIESEN  Coordinates and weights are not yet implemented.
  }

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const
  {
    edgeIds = ArrayView<const gno_t>(Teuchos::null);
    procIds = procIdsConst_(0, numLocalVtx_+1);
    offsets = offsets_(0, numLocalVtx_+1);
    wgts = ArrayView<const scalar_t>(Teuchos::null);

    if (numLocalEdges_){
      if (edgeGnos_.size() == 0)
        edgeIds = edgeGids_(0, numLocalEdges_);
      else
        edgeIds = edgeGnosConst_(0, numLocalEdges_);
    }

    return numLocalEdges_;
  }

  size_t getLocalEdgeList( ArrayView<const lno_t> &edgeIds,
    ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const 
  {
    int nvtx = numLocalVtx_;
    wgts = ArrayView<const scalar_t>(Teuchos::null);

    if (nearEdgeOffsets_.size() > 0){
      if (numNearLocalEdges_ > 0)
        edgeIds = nearEdgeLnos_.view(0, numNearLocalEdges_);
      else
        edgeIds = ArrayView<const lno_t>(Teuchos::null);

      offsets = nearEdgeOffsets_.view(0, nvtx + 1);
    }
    else{
      lno_t *offs = new lno_t [nvtx + 1];
      Z2_LOCAL_MEMORY_ASSERTION(*env_, nvtx+1, offs);
      memset(offs, 0, nvtx+1 * sizeof(lno_t));
      numNearLocalEdges_ = 0;

      if (nvtx > 0){
        for (lno_t i=0; i < nvtx; i++){
          for (lno_t j=offsets_[i]; j < offsets_[i+1]; j++){
            if (procIds_[j] == env_->myRank_){
              offs[i+1]++;
              numNearLocalEdges_++;
            }
          }
        }

        if (numNearLocalEdges_ > 0){
          gno_t *gnos = new gno_t [numNearLocalEdges_];
          Z2_LOCAL_MEMORY_ASSERTION(*env_, numNearLocalEdges_, gnos);
          ArrayView<gno_t> gnoList(gnos, numNearLocalEdges_);

          for (lno_t i=2; i < nvtx; i++){
            offs[i] += offs[i-1];
          }

          for (lno_t i=0; i < nvtx; i++){
            for (lno_t j=offsets_[i]; j < offsets_[i+1]; j++){
              if (procIds_[j] == env_->myRank_){
                gnoList[offs[i]] = edgeGnos_[j];
                offs[i]++;
              }
            }
          }

          lno_t *lnos = new lno_t [numNearLocalEdges_];
          Z2_LOCAL_MEMORY_ASSERTION(*env_, numNearLocalEdges_, lnos);
          ArrayRCP<lno_t> lnoList(lnos, 0, numNearLocalEdges_, true);
          RCP<const idmap_t > idMap = this->getIdentifierMap();

          idMap->lnoTranslate(lnoList.view(0, numNearLocalEdges_), gnoList, 
             TRANSLATE_LIB_TO_APP);

          delete [] gnos;
          nearEdgeLnos_ = arcp_const_cast<const lno_t>(lnoList);
          
          for (lno_t i = nvtx; i > 0; i--){
            offs[i] = offs[i-1];
          }
          offs[0] = 0;
        }
      }

      nearEdgeOffsets_ = arcp(offs, 0, nvtx+1);
    }
    return numNearLocalEdges_;
  }

  ////////////////////////////////////////////////////
  // The Model interface.
  ////////////////////////////////////////////////////

  size_t getLocalNumObjects() const
  {
    return getLocalNumVertices();
  }

  global_size_t getGlobalNumObjects() const
  {
    return getGlobalNumVertices();
  }

  void getGlobalObjectIds(ArrayView<const gno_t> &gnos) const 
  { 
    ArrayView<const scalar_t> xyz, wgts;
    getVertexList(gnos, xyz, wgts);
  }
      }

  int getNumWeights() const { return 0; }

private:

  const MatrixInput<User> *input_;
  const RCP<const Environment > env_;

  ArrayRCP<const gid_t> gids_;
  ArrayRCP<gno_t> gnos_;

  ArrayRCP<const gid_t> edgeGids_;
  ArrayRCP<gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  ArrayRCP<const gno_t> gnosConst_;
  ArrayRCP<const gno_t> edgeGnosConst_;
  ArrayRCP<const int> procIdsConst_;

  size_t numLocalEdges_;
  global_size_t numGlobalEdges_;
  size_t numLocalVtx_;

  bool gidsAreGnos_;

  // For local graphs (graph restricted to local process).  Only
  // create these arrays if required by the algorithm.

  ArrayRCP<const lno_t> nearEdgeLnos_;
  ArrayRCP<const lno_t> nearEdgeOffsets_;
  size_t numNearLocalEdges_;
};

}   // namespace Zoltan2

#endif

