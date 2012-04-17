// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_PartitioningProblem.hpp
    \brief Defines the PartitioningProblem class.
\todo add resetParameters
*/

#ifndef _ZOLTAN2_PARTITIONINGPROBLEM_HPP_
#define _ZOLTAN2_PARTITIONINGPROBLEM_HPP_

#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningAlgorithms.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_GraphModel.hpp>
#include <Zoltan2_IdentifierModel.hpp>

#ifdef HAVE_ZOLTAN2_OVIS
#include <ovis.h>
#endif

namespace Zoltan2{

/*! \brief PartitioningProblem sets up partitioning problems for the user.
 *
 *  The PartitioningProblem is the core of the Zoltan2 partitioning API.
 *  Based on the the user's input and parameters, the PartitioningProblem
 *  sets up a computational Model, and a Solution object.  When the user
 *  calls the solve() method, the PartitioningProblem runs the algorithm,
 *  after which the Solution object may be obtained by the user.
 *  \todo include pointers to examples
 *
 *  The template parameter is the InputAdapter containing the data that
 *  is to be partitioned.
 *
 *  \todo hierarchical partitioning
 *  \todo repartition given an initial solution
 *  \todo follow partitioning with global or local ordering
 *  \todo allow unsetting of part sizes by passing in null pointers
 *  \todo add a parameter by which user tells us there are no self 
 *        edges to be removed.
 */
template<typename Adapter>
class PartitioningProblem : public Problem<Adapter>
{
public:

  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::user_t user_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;

#ifdef HAVE_ZOLTAN2_MPI
  typedef Teuchos::OpaqueWrapper<MPI_Comm> mpiWrapper_t;
#endif

  /*! \brief Destructor
   */
  ~PartitioningProblem() {};

#ifdef HAVE_ZOLTAN2_MPI

  /*! \brief Constructor where MPI communicator can be specified
   */
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p, MPI_Comm comm); 

#endif

  //! \brief Constructor where communicator is the Teuchos default.
  PartitioningProblem(Adapter *A, Teuchos::ParameterList *p) ;


  //!  \brief Reset the parameter list.

  void resetParameterList(Teuchos::ParameterList *p);

  //!  \brief Direct the problem to create a solution.
  //
  //    \param updateInputData   If true this indicates that either
  //          this is the first attempt at solution, or that we
  //          are computing a new solution and the input data has
  //          changed since the previous solution was computed.  
  //          By input data we mean coordinates, topology, or weights.
  //          If false, this indicates that we are computing a
  //          new solution using the same input data was used for
  //          the previous solution, even though the parameters
  //          may have been changed.
  //
  //  For the sake of performance, we ask the caller to set \c updateInputData
  //  to false if he/she is computing a new solution using the same input data,
  //  but different problem parameters, than that which was used to compute 
  //  the most recent solution.

  void solve(bool updateInputData=true );
 
  //!  \brief Get the solution to the problem.
  //
  //   \return  a reference to the solution to the most recent solve().

  const PartitioningSolution<Adapter> &getSolution() {
    return *(solution_.getRawPtr());
  };

  /*! \brief Set or reset relative sizes for the partitions that Zoltan2 will create.
   *
   *  \param len  The size of the \c partIds and \c partSizes lists
   *  \param partIds   A list of \c len partition identifiers.  Partition
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of \c len relative sizes corresponding to
   *           the \c partIds.
   *  \param makeCopy  If true, Zoltan2 will make a copy of the ids and sizes
   *      that are provided in this call.  If false, Zoltan2 will just save
   *      the pointers to to the caller's lists.  If the pointers will remain
   *      remain valid throughout the lifetime of the PartitioningProblem,
   *      and memory use is an issue, then set makeCopy to false.  By default,
   *      Zoltan2 will copy the caller's list of ids and sizes.
   *
   * A given partid should only be provided once across the application.
   * Duplicate partIds will generate a std::runtime_error exception when
   * the PartitioningSolution is created.  Part
   * ids that are omitted will be assigned the average of the sizes that
   * have been specified.
   *
   * Subsequent calls to setPartSizes will replace the list of part ids
   * and part sizes provided previously.
   * 
   * If the application has set multiple weights per object, then the
   * part sizes supplied in this method are applied to the first weight.
   *
   * Zoltan2 assumes that uniform partition sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCritiera.
   *
   * \todo A user should be able to give us one set of part sizes
   *            that applies to all weight dimensions.  Right now
   *            for each weight dimension that does not have
   *            uniform part sizes, the user has to give us the
   *            part sizes once for each.
   */

  void setPartSizes(int len, partId_t *partIds, scalar_t *partSizes, 
    bool makeCopy=true) 
  { 
    setPartSizesForCritiera(0, len, partIds, partSizes, makeCopy);
  }

  /*! \brief Set or reset the relative sizes (per weight) for the partitions 
   *    that Zoltan2 will create.
   *
   *  \param criteria the criteria (weight dimension) for which these 
   *      part sizes apply.  Criteria range from zero to one less than
   *     the number of weights per object specified in the 
   *     caller's InputAdapter.
   *  \param len  The size of the \c partIds and \c partSizes lists
   *  \param partIds   A list of \c len partition identifiers.  Partition
   *           identifiers range from zero to one less than the global
   *           number of identifiers.  
   *  \param partSizes  A list of \c len relative sizes corresponding to
   *           the \c partIds.
   *  \param makeCopy  If true, Zoltan2 will make a copy of the ids and sizes
   *      that are provided in this call.  If false, Zoltan2 will just save
   *      the pointers to to the caller's lists.  If the pointers will remain
   *      remain valid throughout the lifetime of the PartitioningProblem,
   *      and memory use is an issue, then set makeCopy to false.  By default,
   *      Zoltan2 will copy the caller's list of ids and sizes.
   *
   * A given partid should only be provided once across the application.
   * Duplicate partIds will generate a std::runtime_error exception when
   * the PartitioningSolution is created.  Part
   * ids that are omitted will be assigned the average of the sizes that
   * have been specified.
   *
   * Subsequent calls to setPartSizes for the same criteria will replace 
   * the list of part ids and part sizes provided for that criteria previously.
   *
   * Zoltan2 assumes that uniform partition sizes are desired by the caller,
   * unless specified otherwise in a call to setPartSizes or 
   * setPartSizesForCritiera.
   */

  void setPartSizesForCritiera(int criteria, int len, partId_t *partIds, 
    scalar_t *partSizes, bool makeCopy=true) ;

private:
  void initializeProblem();

  void createPartitioningProblem(bool newData);

  RCP<PartitioningSolution<Adapter> > solution_;

  RCP<Comm<int> > problemComm_;
  RCP<const Comm<int> > problemCommConst_;

#ifdef HAVE_ZOLTAN2_MPI
  MPI_Comm mpiComm_;
#endif

  InputAdapterType inputType_;
  ModelType modelType_;
  modelFlag_t graphFlags_;
  modelFlag_t idFlags_;
  modelFlag_t coordFlags_;
  std::string algorithm_;

  int numberOfWeights_;

  // Suppose Array<partId_t> partIds = partIds_[w].  If partIds.size() > 0
  // then the user supplied part sizes for weight index "w", and the sizes
  // corresponding to the Ids in partIds are partSizes[w].
  //
  // If numberOfWeights_ >= 0, then there is an Id and Sizes array for
  // for each weight.  Otherwise the user did not supply object weights,
  // but they can still specify part sizes. 
  // So numberOfCriteria_ is numberOfWeights_ or one, whichever is greater.

  ArrayRCP<ArrayRCP<partId_t> > partIds_;
  ArrayRCP<ArrayRCP<scalar_t> > partSizes_;
  int numberOfCriteria_;

  // Number of parts to be computed at each level in hierarchical partitioning.
  
  ArrayRCP<int> levelNumberParts_;
  bool hierarchical_;
};
////////////////////////////////////////////////////////////////////////

#ifdef HAVE_ZOLTAN2_MPI
template <typename Adapter>
  PartitioningProblem<Adapter>::PartitioningProblem(Adapter *A, 
    ParameterList *p, MPI_Comm comm):
      Problem<Adapter>(A,p,comm), solution_(),
      problemComm_(), problemCommConst_(),
      inputType_(InvalidAdapterType), modelType_(InvalidModel), 
      graphFlags_(), idFlags_(), coordFlags_(), algorithm_(),
      numberOfWeights_(), partIds_(), partSizes_(), 
      numberOfCriteria_(), levelNumberParts_(), hierarchical_(false)
{
  initializeProblem();
}
#endif

template <typename Adapter>
  PartitioningProblem<Adapter>::PartitioningProblem(Adapter *A, 
    ParameterList *p):
      Problem<Adapter>(A,p), solution_(),
      problemComm_(), problemCommConst_(),
      inputType_(InvalidAdapterType), modelType_(InvalidModel), 
      graphFlags_(), idFlags_(), coordFlags_(), algorithm_(),
      numberOfWeights_(), 
      partIds_(), partSizes_(), numberOfCriteria_(), 
      levelNumberParts_(), hierarchical_(false)
{
  initializeProblem();
}

template <typename Adapter>
  void PartitioningProblem<Adapter>::initializeProblem()
{
  HELLO;

#ifdef HAVE_ZOLTAN2_OVIS
  ovis_enabled(this->comm_->getRank());
#endif

  // Create a copy of the user's communicator.  

  problemComm_ = this->comm_->duplicate();
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);

#ifdef HAVE_ZOLTAN2_MPI

  // TPLs may want an MPI communicator

  Comm<int> *c = problemComm_.getRawPtr();
  Teuchos::MpiComm<int> *mc = dynamic_cast<Teuchos::MpiComm<int> *>(c);
  if (mc){
    RCP<const mpiWrapper_t> wrappedComm = mc->getRawMpiComm();
    mpiComm_ = (*wrappedComm.getRawPtr())();
  }
  else{
    mpiComm_ = MPI_COMM_SELF;   // or would this be an error?
  }

#endif

  // Number of criteria is number of user supplied weights if non-zero.
  // Otherwise it is 1 and uniform weight is implied.

  numberOfWeights_ = this->inputAdapter_->getNumberOfWeightsPerObject();

  numberOfCriteria_ = (numberOfWeights_ > 1) ? numberOfWeights_ : 1;

  inputType_ = this->inputAdapter_->inputAdapterType();

  // The Caller can specify part sizes in setPartSizes().  If he/she
  // does not, the part size arrays are empty.

  ArrayRCP<partId_t> *noIds = new ArrayRCP<partId_t> [numberOfCriteria_];
  ArrayRCP<scalar_t> *noSizes = new ArrayRCP<scalar_t> [numberOfCriteria_];

  partIds_ = arcp(noIds, 0, numberOfCriteria_, true);
  partSizes_ = arcp(noSizes, 0, numberOfCriteria_, true);
}

template <typename Adapter>
  void PartitioningProblem<Adapter>::setPartSizesForCritiera(
    int criteria, int len, partId_t *partIds, scalar_t *partSizes, bool makeCopy) 
{
  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid length", 
    len>= 0, BASIC_ASSERTION);

  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid criteria", 
    criteria >= 0 && criteria < numberOfWeights_, BASIC_ASSERTION);

  if (len == 0){
    partIds_[criteria] = ArrayRCP<partId_t>();
    partSizes_[criteria] = ArrayRCP<scalar_t>();
    return;
  }

  this->env_->localInputAssertion(__FILE__, __LINE__, "invalid arrays", 
    partIds && partSizes, BASIC_ASSERTION);

  // The global validity of the partIds and partSizes arrays is performed
  // by the PartitioningSolution, which computes global part distribution and
  // part sizes.

  partId_t *z2_partIds = partIds;
  scalar_t *z2_partSizes = partSizes;
  bool own_memory = false;

  if (makeCopy){
    z2_partIds = NULL;
    z2_partIds = new partId_t [len];
    this->env_->localMemoryAssertion(__FILE__, __LINE__, len, z2_partIds);
    z2_partSizes = NULL;
    z2_partSizes = new scalar_t [len];
    this->env_->localMemoryAssertion(__FILE__, __LINE__, len, z2_partSizes);
  }

  partIds_[criteria] = arcp(z2_partIds, 0, len, own_memory);
  partSizes_[criteria] = arcp(z2_partSizes, 0, len, own_memory);
}

template <typename Adapter>
void PartitioningProblem<Adapter>::solve(bool updateInputData)
{
  HELLO;

  // Create the computational model.

  createPartitioningProblem(updateInputData);

  // TODO: If hierarchical_

  // Create the solution. The algorithm will query the Solution
  //   for part and weight information. The algorithm will
  //   update the solution with part assignments and quality
  //   metrics.  The Solution object itself will convert our internal
  //   global numbers back to application global Ids if needed.

  RCP<const IdentifierMap<user_t> > idMap = 
    this->baseModel_->getIdentifierMap();

  PartitioningSolution<Adapter> *soln = NULL;

  try{
    soln = new PartitioningSolution<Adapter>( 
      this->envConst_, problemCommConst_, idMap, numberOfWeights_, 
      partIds_.view(0, numberOfCriteria_), 
      partSizes_.view(0, numberOfCriteria_));
  }
  Z2_FORWARD_EXCEPTIONS;

  solution_ = rcp(soln);

  // Call the algorithm

  try {
    if (algorithm_ == string("scotch")){
      AlgPTScotch<Adapter>(this->envConst_, problemComm_,
#ifdef HAVE_ZOLTAN2_MPI
        mpiComm_,
#endif
        this->graphModel_, solution_);
    }
    else if (algorithm_ == string("block")){
      AlgBlock<Adapter>(this->envConst_, problemComm_,
        this->identifierModel_, solution_);
    }
    else if (algorithm_ == string("rcb")){
      AlgRCB<Adapter>(this->envConst_, problemComm_,
        this->coordinateModel_, solution_);
    }
    else{
      throw std::logic_error("partitioning algorithm not supported yet");
    }
  }
  Z2_FORWARD_EXCEPTIONS;

#ifdef HAVE_ZOLTAN2_MPI

  // The algorithm may have changed the communicator.  Change it back.

  RCP<const mpiWrapper_t > wrappedComm = rcp(new mpiWrapper_t(mpiComm_));
  problemComm_ = rcp(new Teuchos::MpiComm<int>(wrappedComm));
  problemCommConst_ = rcp_const_cast<const Comm<int> > (problemComm_);

#endif

}

template <typename Adapter>
void PartitioningProblem<Adapter>::createPartitioningProblem(bool newData)
{
  HELLO;
  using std::string;
  using Teuchos::ParameterList;

  // A Problem object may be reused.  The input data may have changed and
  // new parameters or part sizes may have been set.
  //
  // Save these values in order to determine if we need to create a new model.

  ModelType previousModel = modelType_;
  modelFlag_t previousGraphModelFlags = graphFlags_;
  modelFlag_t previousIdentifierModelFlags = idFlags_;
  modelFlag_t previousCoordinateModelFlags = coordFlags_;

  modelType_ = InvalidModel;
  graphFlags_.reset();
  idFlags_.reset();
  coordFlags_.reset();

  ////////////////////////////////////////////////////////////////////////////
  // It's possible at this point that the Problem may want to
  // add problem parameters to the parameter list in the Environment. 
  //
  // Since the parameters in the Environment have already been
  // validated in its constructor, a new Environment must be created:
  ////////////////////////////////////////////////////////////////////////////
  // Teuchos::RCP<const Teuchos::Comm<int> > oldComm = this->env_->comm_;
  // const ParameterList &oldParams = this->env_->getUnvalidatedParameters();
  // 
  // ParameterList newParams = oldParams;
  // newParams.set("new_parameter", "new_value");
  // 
  // ParameterList &newPartParams = newParams.sublist("partitioning");
  // newPartParams.set("new_partitioning_parameter", "its_value");
  // 
  // this->env_ = rcp(new Environment(newParams, oldComm));
  ////////////////////////////////////////////////////////////////////////////

  Environment &env = *(this->env_);
  ParameterList &pl = env.getParametersNonConst();

  bool isSet;
  string defString("default");

  // Did the user specify a computational model?

  string model(defString);
  getParameterValue(pl, "partitioning", "model", isSet, model);

  // Did the user specify an algorithm?

  string algorithm(defString);
  getParameterValue(pl, "partitioning", "algorithm", isSet, algorithm);

  // Possible algorithm requirements that must be conveyed to the model:

  bool needConsecutiveGlobalIds = false;
  bool removeSelfEdges= false;

  ///////////////////////////////////////////////////////////////////
  // Determine algorithm, model, and algorithm requirements.  This
  // is a first pass.  Feel free to change this and add to it.
  
  if (algorithm != defString){
    // Figure out the model required by the algorithm
    if (algorithm == string("block") ||
        algorithm == string("random") ||
        algorithm == string("cyclic") ){

      modelType_ = IdentifierModelType;
      algorithm_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == string("rcb") ||
             algorithm == string("rib") ||
             algorithm == string("hsfc")){

      modelType_ = CoordinateModelType;
      algorithm_ = algorithm;
    }
    else if (algorithm == string("metis") ||
             algorithm == string("parmetis") ||
             algorithm == string("scotch") ||
             algorithm == string("ptscotch")){

      modelType_ = GraphModelType;
      algorithm_ = algorithm;
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
    }
    else if (algorithm == string("patoh") ||
             algorithm == string("phg")){

      if ((modelType_ != GraphModelType) &&
          (modelType_ != HypergraphModelType) ){
        modelType_ = HypergraphModelType;
      }
      algorithm_ = algorithm;
      needConsecutiveGlobalIds = true;
    }
    else{
      // Parameter list should ensure this does not happen.
      throw std::logic_error("parameter list algorithm is invalid");
    }
  }
  else if (model != defString){
    // Figure out the algorithm suggested by the model.
    if (model == string("hypergraph")){
      modelType_ = HypergraphModelType;
      if (problemComm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
      needConsecutiveGlobalIds = true;
    }
    else if (model == string("graph")){
      modelType_ = GraphModelType;
#ifdef HAVE_ZOLTAN2_SCOTCH
      if (problemComm_->getSize() > 1)
        algorithm_ = string("ptscotch"); 
      else
        algorithm_ = string("scotch"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
#ifdef HAVE_ZOLTAN2_PARMETIS
      if (problemComm_->getSize() > 1)
        algorithm_ = string("parmetis"); 
      else
        algorithm_ = string("metis"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#else
      if (problemComm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
      removeSelfEdges = true;
      needConsecutiveGlobalIds = true;
#endif
#endif
    }
    else if (model == string("geometry")){
      modelType_ = CoordinateModelType;
      algorithm_ = string("rcb");
    }
    else if (model == string("ids")){
      modelType_ = IdentifierModelType;
      algorithm_ = string("block");
      needConsecutiveGlobalIds = true;
    }
    else{
      // Parameter list should ensure this does not happen.
      env.localBugAssertion(__FILE__, __LINE__, 
        "parameter list model type is invalid", 1, BASIC_ASSERTION);
    }
  }
  else{   
    // Determine an algorithm and model suggested by the input type.
    //   TODO: this is a good time to use the time vs. quality parameter
    //     in choosing an algorithm, and setting some parameters

    if (inputType_ == MatrixAdapterType){
      modelType_ = HypergraphModelType;
      if (problemComm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
    }
    else if (inputType_ == GraphAdapterType ||
        inputType_ == MeshAdapterType){
      modelType_ = GraphModelType;
      if (problemComm_->getSize() > 1)
        algorithm_ = string("phg"); 
      else
        algorithm_ = string("patoh"); 
    }
    else if (inputType_ == CoordinateAdapterType){
      modelType_ = CoordinateModelType;
      algorithm_ = string("rcb");
    }
    else if (inputType_ == VectorAdapterType ||
             inputType_ == IdentifierAdapterType){
      modelType_ = IdentifierModelType;
      algorithm_ = string("block");
    }
    else{
      // This should never happen
      throw std::logic_error("input type is invalid");
    }
  }

  // Hierarchical partitioning?

  Array<int> valueList;
  getParameterValue(pl, "partitioning", "topology", isSet, valueList);

  if (isSet){
    if (!Zoltan2::noValuesAreInRangeList(valueList)){
      int *n = new int [valueList.size() + 1];
      levelNumberParts_ = arcp(n, 0, valueList.size() + 1, true);
      int procsPerNode = 1;
      for (int i=0; i < valueList.size(); i++){
        levelNumberParts_[i+1] = valueList[i];
        procsPerNode *= valueList[i];
      }
      // Number of parts in the first level
      levelNumberParts_[0] = env.numProcs_ / procsPerNode;

      if (env.numProcs_ % procsPerNode > 0)
        levelNumberParts_[0]++;
    }
  }
  else{
    levelNumberParts_.clear();
  }

  hierarchical_ = levelNumberParts_.size() > 0;

  // Object to be partitioned? (rows, columns, etc)

  string objectOfInterest(defString);
  getParameterValue(pl, "partitioning", "objects", isSet, objectOfInterest);

  ///////////////////////////////////////////////////////////////////
  // Set model creation flags, if any.

  if (modelType_ == GraphModelType){

    // Any parameters in the graph sublist?

    string symParameter(defString);
    getParameterValue(pl, "partitioning", "graph", "symmetrize_input",
      isSet, symParameter);

    if (isSet){
      if (symParameter == string("transpose"))
        graphFlags_.set(SYMMETRIZE_INPUT_TRANSPOSE);
      else if (symParameter == string("bipartite"))
        graphFlags_.set(SYMMETRIZE_INPUT_BIPARTITE);
    } 

    int sgParameter = 0;
    getParameterValue(pl, "partitioning", "graph", "subset_graph",
      isSet, sgParameter);

    if (isSet && (sgParameter == 1))
        graphFlags_.set(GRAPH_IS_A_SUBSET_GRAPH);

    // Any special behaviors required by the algorithm?
    
    if (removeSelfEdges)
      graphFlags_.set(SELF_EDGES_MUST_BE_REMOVED);

    if (needConsecutiveGlobalIds)
      graphFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);

    // How does user input map to vertices and edges?

    if (inputType_ == MatrixAdapterType){
      if (objectOfInterest == defString ||
          objectOfInterest == string("matrix_rows") )
        graphFlags_.set(VERTICES_ARE_MATRIX_ROWS);
      else if (objectOfInterest == string("matrix_columns"))
        graphFlags_.set(VERTICES_ARE_MATRIX_COLUMNS);
      else if (objectOfInterest == string("matrix_nonzeros"))
        graphFlags_.set(VERTICES_ARE_MATRIX_NONZEROS);
    }

    else if (inputType_ == MeshAdapterType){
      if (objectOfInterest == defString ||
          objectOfInterest == string("mesh_nodes") )
        graphFlags_.set(VERTICES_ARE_MESH_NODES);
      else if (objectOfInterest == string("mesh_elements"))
        graphFlags_.set(VERTICES_ARE_MESH_ELEMENTS);
    } 
  }
  else if (modelType_ == IdentifierModelType){

    // Any special behaviors required by the algorithm?
    
    if (needConsecutiveGlobalIds)
      idFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  }
  else if (modelType_ == CoordinateModelType){

    // Any special behaviors required by the algorithm?
    
    if (needConsecutiveGlobalIds)
      coordFlags_.set(IDS_MUST_BE_GLOBALLY_CONSECUTIVE);
  }


  if (  newData ||
       (modelType_ != previousModel) ||
       (graphFlags_ != previousGraphModelFlags) ||
       (coordFlags_ != previousCoordinateModelFlags) ||
       (idFlags_ != previousIdentifierModelFlags) ) {

    // Create the computational model.
    // Models are instantiated for base input adapter types (mesh,
    // matrix, graph, and so on).  We pass a pointer to the input
    // adapter, cast as the base input type.

    typedef typename Adapter::base_adapter_t base_adapter_t;

    switch (modelType_) {

    case GraphModelType:
      this->graphModel_ = rcp(new GraphModel<base_adapter_t>(
        this->baseInputAdapter_, this->envConst_, problemComm_, graphFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
        this->graphModel_);

      break;

    case HypergraphModelType:
      break;
  
    case CoordinateModelType:
      this->coordinateModel_ = rcp(new CoordinateModel<base_adapter_t>(
        this->baseInputAdapter_, this->envConst_, problemComm_, coordFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
        this->coordinateModel_);
      break;

    case IdentifierModelType:
      this->identifierModel_ = rcp(new IdentifierModel<base_adapter_t>(
        this->baseInputAdapter_, this->envConst_, problemComm_, idFlags_));

      this->baseModel_ = rcp_implicit_cast<const Model<base_adapter_t> >(
        this->identifierModel_);
      break;

    default:
      cout << __func__ << " Invalid model" << modelType_ << endl;
      break;
    }
  }
}

}  // namespace Zoltan2
#endif
