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

#include <vector>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Hashtable.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Zoltan2_Model.hpp>
#include <Zoltan2_XpetraCrsMatrixInput.hpp>
#include <Zoltan2_IdentifierMap.hpp>
#include <Xpetra_Map.hpp>

namespace Zoltan2 {

/*! Zoltan2::GraphModel
    \brief GraphModel defines the interface required for graph models.  

    The constructor of the GraphModel can be a global call, requiring
    all processes in the application to call it.  The rest of the
    method should be local methods.

    The template parameter is an Input Adapter.  Input adapters are
    templated on the basic user input type.

    TODO: GraphModel parameters: 
                    consecutive IDs required
                    base required
                    remove self edges
*/
template <typename Adapter>
class GraphModel : public Model<Adapter>
{
private:

public:

  typedef typename Adapter::scalar_t  scalar_t;
  typedef typename Adapter::gno_t     gno_t;
  typedef typename Adapter::lno_t     lno_t;
  typedef typename Adapter::gid_t     gid_t;
  typedef typename Adapter::lid_t     lid_t;
  typedef typename Adapter::node_t    node_t;
  typedef typename Adapter::user_t    user_t;
  
  GraphModel(){
    //std::cout <<"GRAPH MODEL base" << std::endl;
  }

  /*! Returns the number vertices on this process.
   */
  size_t getLocalNumVertices() const { return 0; }

  /*! Returns the global number vertices.
   */
  global_size_t getGlobalNumVertices() const { return 0; }

  /*! Returns the number edges on this process.
   */
  size_t getLocalNumEdges() const { return 0; }

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

  /*! Returns the base ID, typically 0 or 1.
   */
  gno_t getIndexBase() const {return 0;}

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
   */

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const {
      return 0; }

  /*! Sets pointers to this process' edge (neighbor) global Ids.
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

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds,
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const { return 0; }

  /*! Obtain a view of the edge Ids of the input vertex.
      \param Id  is the global Id for a vertex on this process.
      \param edgeIds on return will point to the list of edge neighbors.
      \param procIds on return holds the list of each process owning the
        corresponding neighbor vertex. 
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.
  
      This method is defined for convenience when obtaining the
      neighbors of a vertex.  It is not efficient to call this method
      many times in a loop, due to the construction and destruction of
      ArrayViews.  Call getEdgeList instead.
   */
  size_t getVertexGlobalEdge( gno_t Id, 
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const scalar_t> &wgts) const { return 0; }
   
  /*! Obtain a view of the edge Ids of the input vertex.
      \param localRef  is the local id associated with vertex.  Local ids
        are consecutive, begin at 0 and follow the order of vertices returned
        by getVertexList.
      \param edgeIds on return will point to the list of edge neighbor global
         Ids.
      \param procIds on return holds the list of each process owning the
        corresponding neighbor vertex.
      \param wgts on return points to the weights, if any, associated with the
         edges. Weights are listed by edge by weight component.
      \return The number of ids in the edgeId list.

      This method is defined for convenience when obtaining the 
      neighbors of a vertex.  It is not efficient to call this method 
      many times in a loop, due to the construction and destruction of
      ArrayViews.  Call getEdgeList instead.
   */
  size_t getVertexLocalEdge( lno_t localRef,
    ArrayView<const gno_t> &edgeIds, ArrayView<const int> &procIds,
    ArrayView<const scalar_t> &wgts) const { return 0; }
};

////////////////////////////////////////////////////////////////
// Graph model derived from XpetraCrsMatrixInput.
//    We know that Xpetra input does not need an IdentifierMap
//    unless consecutive global Ids are required and the Xpetra
//    input is not consecutive.
//
//   TODO Xpetra matrices need not have consecutive IDs, so
//               add use if IdentifierMap to ensure IDs are
//               consecutive
////////////////////////////////////////////////////////////////

/*! Zoltan2::GraphModel<XpetraCrsMatrixInput>
    \brief A (partial) specialization of GraphModel
           for a Zoltan2::XpetraCrsMatrixInput object.
*/

template <typename User>
class GraphModel<XpetraCrsMatrixInput<User> >
{
public:

  typedef typename XpetraCrsMatrixInput<User>::scalar_t  scalar_t;
  typedef typename XpetraCrsMatrixInput<User>::gno_t     gno_t;
  typedef typename XpetraCrsMatrixInput<User>::lno_t     lno_t;
  typedef typename XpetraCrsMatrixInput<User>::gid_t     gid_t;
  typedef typename XpetraCrsMatrixInput<User>::lid_t     lid_t;
  typedef typename XpetraCrsMatrixInput<User>::node_t    node_t;

  /*! Constructor
   *  All processes in the communicator must call the constructor.
   *  \param  inputAdapter  an encapsulation of the user data
   *  \param  comm          the problem communicator
   *  \param  env           environment (library configuration settings)
   *  \param  env           environment (library configuration settings)
   *  \param  consecutiveIdsRequired  set to true if the algorithm or
   *           third party library requires consecutive global vertex Ids.
   */
  GraphModel(
    const RCP<const XpetraCrsMatrixInput<User> > &inputAdapter,
    const RCP<const Comm<int> > &comm, const RCP<const Environment> &env,
    bool consecutiveIdsRequired=false, bool removeSelfEdges=false) :
      input_(inputAdapter), rowMap_(inputAdapter->getMatrix()->getRowMap()),
      colMap_(inputAdapter->getMatrix()->getColMap()), comm_(comm), env_(env),
      gnos_(), edgeGnos_(), procIds_(), offsets_(),
      numLocalEdges_(), numGlobalEdges_(0), numLocalVtx_()
  {
    gno_t const *vtxIds=NULL, *nborIds=NULL;
    lno_t const  *offsets=NULL, *lids=NULL; 
    try{
      numLocalVtx_ = input_->getRowListView(vtxIds, lids, offsets, nborIds);
    }
    catch (std::exception &e)
      Z2_THROW_ZOLTAN2_ERROR(env_, e);

    numLocalEdges_ = offsets[numLocalVtx_];
    size_t numOffsets = numLocalVtx_ + 1;

    ArrayView<gno_t> av1(const_cast<gno_t *>(vtxIds), numLocalVtx_);
    gnos_ = av1.getConst();  // to make ArrayView<const gno_t>

    ArrayRCP<lno_t> tmpOffsets;
    ArrayRCP<gno_t> tmpEdges;
    lno_t nSelfEdges = 0;

    if (removeSelfEdges) {

      lno_t *offArray = new lno_t [numOffsets];
      Z2_LOCAL_MEMORY_ASSERTION(*env, numOffsets, offArray);
      gno_t *edArray = new gno_t [numLocalEdges_];
      Z2_LOCAL_MEMORY_ASSERTION(*env, numLocalEdges_, !numLocalEdges_||edArray);

      for (lno_t i=0; i < numLocalVtx_; i++){

        offArray[i] = offsets[i] - nSelfEdges;

        for (lno_t j = offsets[i]; j < offsets[i+1]; j++) {
          if (gnos_[i] == nborIds[j]) { // self edge; remove it
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
        tmpOffsets = arcp(offArray, 0, numLocalVtx_+1, true);
        tmpEdges = arcp(edArray, 0, numLocalEdges_, true);
      }
      else{
        delete [] offArray;
        if (numLocalEdges_) delete [] edArray;
      }
    }
  
    if (nSelfEdges == 0){
      offsets_ = arcp(const_cast<lno_t *>(offsets), 0, numOffsets, false);
      edgeGnos_ = arcp(const_cast<gno_t *>(nborIds), 0, numLocalEdges_, false);
    }
    else{
      offsets_ = tmpOffsets;
      edgeGnos_ =  tmpEdges;
    }

    Teuchos::reduceAll<int, size_t>(*comm, Teuchos::REDUCE_SUM, 1,
      &numLocalEdges_, &numGlobalEdges_);

    RCP<Array<int> > procBuf =  rcp(new Array<int>(numLocalEdges_));
    procIds_ = arcp(procBuf);

    try{
      ArrayView<const gno_t> gnoView = 
        edgeGnos_.view(0,numLocalEdges_).getConst();
      rowMap_->getRemoteIndexList(gnoView, procIds_.view(0, numLocalEdges_));
    }
    catch (std::exception &e){
      Z2_THROW_ZOLTAN2_ERROR(env_, e);
    }
  }

  //!  Destructor
  ~GraphModel() {
  }

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

  gno_t getIndexBase() const
  {
  }

  size_t getVertexList( ArrayView<const gno_t> &Ids,
    ArrayView<const scalar_t> &xyz, ArrayView<const scalar_t> &wgts) const
  {
    size_t n = gnos_.size();
    Ids = gnos_.view(0, n);
    return n;
  }

  size_t getEdgeList( ArrayView<const gno_t> &edgeIds, 
    ArrayView<const int> &procIds, ArrayView<const lno_t> &offsets,
    ArrayView<const scalar_t> &wgts) const
  {
    edgeIds = edgeGnos_.view(0,numLocalEdges_).getConst();
    procIds = procIds_.view(0, numLocalEdges_).getConst();
    offsets = offsets_.view(0, numLocalVtx_+1).getConst();

    return numLocalEdges_;
  }

  size_t getVertexGlobalEdge( gno_t Id, ArrayView<const gno_t> &edgeId,
    ArrayView<const int> &procId, ArrayView<const scalar_t> &wgts) const
  {
    if (rowMap_->isNodeGlobalElement(Id)){
      return getVertexLocalEdge(rowMap_->getLocalElement(Id), edgeId, procId, wgts);
    }
    else{
      throw std::runtime_error("Global Id not on this process");
    }
  }

  size_t getVertexLocalEdge( lno_t lno, ArrayView<const gno_t> &edgeId,
    ArrayView<const int> &procId, ArrayView<const scalar_t> &wgts) const
  { 
    Z2_LOCAL_INPUT_ASSERTION(*env_, "invalid local id",
      lno >= 0 && lno < gnos_.size(), BASIC_ASSERTION);

    lno_t thisVtx =  offsets_[lno];
    lno_t nextVtx = (lno < gnos_.size()-1) ? offsets_[lno+1] : numLocalEdges_;
    size_t nEdges = nextVtx - thisVtx;

    edgeId = edgeGnos_.view(thisVtx, nEdges).getConst();
    procId = procIds_.view(thisVtx, nEdges).getConst();
    return nEdges;
  }

private:

  RCP<const XpetraCrsMatrixInput<User> > input_;
  const RCP<const Xpetra::Map<lno_t, gno_t> > rowMap_;
  const RCP<const Xpetra::Map<lno_t, gno_t> > colMap_;
  RCP<const Teuchos::Comm<int> > comm_;
  RCP<const Environment > env_;

  ArrayView<const gno_t> gnos_;
  ArrayRCP<const gno_t> edgeGnos_;
  ArrayRCP<int> procIds_;
  ArrayRCP<const lno_t> offsets_;

  // Transpose is required only if vertices are columns.
  // KDDKDD ??  We won't form an actual transpose, will we?
  // KDDKDD RCP<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> > inputTranspose;

  global_size_t numLocalEdges_;
  global_size_t numGlobalEdges_;
  size_t numLocalVtx_;
};

}   // namespace Zoltan2

#endif
