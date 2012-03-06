// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraVectorInput.hpp
    \brief Defines the XpetraVectorInput adapter class.
*/

#ifndef _ZOLTAN2_XPETRAVECTORINPUT_HPP_
#define _ZOLTAN2_XPETRAVECTORINPUT_HPP_

#include <Zoltan2_VectorInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_StridedInput.hpp>
#include <Zoltan2_Util.hpp>

#include <Xpetra_EpetraVector.hpp>
#include <Xpetra_TpetraVector.hpp>

namespace Zoltan2 {


/*!  \brief Provides access for Zoltan2 to an Xpetra::Vector.

    The template parameter is the user's input data type, which can be:
   \li Epetra_Vector
   \li Tpetra::Vector
   \li Xpetra::Vector
*/

template <typename User>
class XpetraVectorInput : public VectorInput<User> {
public:

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef VectorInput<User>       base_adapter_t;
  typedef User user_t;

  typedef Xpetra::Vector<
    scalar_t, lno_t, gno_t, node_t> x_vector_t;
  typedef Xpetra::TpetraVector<
    scalar_t, lno_t, gno_t, node_t> xt_vector_t;
  typedef Xpetra::EpetraVector xe_vector_t;
#endif

  /*! \brief Destructor
   */
  ~XpetraVectorInput() { }

  /*! \brief Constructor   
   *
   *  \param invector  the user's Xpetra, Tpetra or Epetra Vector object
   *  \param numWeights the number of weights per element, which may be zero
   *                or greater
   *  \param weights  numWeights pointers to arrays of weights
   *  \param weightStrides  a list of numWeights strides for the weights
   *        arrays. The n'th weight for multivector element k is to be found
   *               at weights[n][k*weightStrides[n]].  If weightStrides
   *              is NULL, it is assumed all strides are one.
   *
   *  The values pointed to the arguments must remain valid for the
   *  lifetime of this InputAdapter.
   */
  XpetraVectorInput( const RCP<const User> &invector, int numWeights, 
    const scalar_t * * const weights, int *weightStrides);

  /*! \brief Access to the xpetra-wrapped vector
   */

  const RCP<const x_vector_t> &getVector() const
  {
    return vector_;
  }

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  std::string inputAdapterName()const { return std::string("XpetraVector");}

  size_t getLocalNumberOfObjects() const { return getLocalLength();}

  int getNumberOfWeightsPerObject() const { return numWeights_;}

  ////////////////////////////////////////////////////
  // The VectorInput interface.
  ////////////////////////////////////////////////////

  int getNumberOfVectors() const { return 1; }

  int getNumberOfWeights() const {return numWeights_;}

  size_t getLocalLength() const {return vector_->getLocalLength();}
  
  size_t getGlobalLength() const {return vector_->getGlobalLength();}

  size_t getVector(const gid_t *&Ids, const scalar_t *&elements, 
    int &stride) const;

  size_t getVector(int vectorNumber, const gid_t *&Ids, 
    const scalar_t *&elements, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid vector",
      vectorNumber==0, BASIC_ASSERTION);

    return getVector(Ids, elements, stride);
  }

  size_t getVectorWeights(int dim, const scalar_t *&weights, int &stride) const
  {
    env_->localInputAssertion(__FILE__, __LINE__, "invalid dimension",
      dim >= 0 && dim < numWeights_, BASIC_ASSERTION);

    size_t length;

    weights_[dim]->getStridedList(length, weights, stride);

    return length;
  }

  template <typename User2>
    size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<User2> &solution) const;

private:

  RCP<const User> invector_;
  RCP<const x_vector_t> vector_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > map_;
  RCP<Environment> env_;
  lno_t base_;

  int numWeights_;
  Array<RCP<StridedInput<lno_t, scalar_t> > > weights_;
};

////////////////////////////////////////////////////////////////
// Definitions
////////////////////////////////////////////////////////////////
  
template <typename User>
  XpetraVectorInput<User>::XpetraVectorInput(const RCP<const User> &invector, 
    int numWeights, const scalar_t * * const weights, int *weightStrides):
      invector_(invector), vector_(), map_(),
      env_(rcp(new Environment)), base_(),
      numWeights_(numWeights), weights_(numWeights)
{
  typedef StridedInput<lno_t, scalar_t> input_t;

  vector_ = XpetraTraits<User>::convertToXpetra(invector);
  map_ = vector_->getMap();
  base_ = map_->getIndexBase();

  size_t length = vector_->getLocalLength();

  if (length > 0 && numWeights > 0){
    int stride = 1;
    for (int w=0; w < numWeights; w++){
      if (weightStrides)
        stride = weightStrides[w];
      weights_[w] = rcp<input_t>(new input_t(
        ArrayView<const scalar_t>(weights[w], stride*length), stride));
    }
  }
}

template <typename User>
  size_t XpetraVectorInput<User>::getVector(const gid_t *&Ids, 
    const scalar_t *&elements, int &stride) const
{
  stride = 1;
  elements = NULL;
  const x_vector_t *vec =  vector_.get();

  if (map_->lib() == Xpetra::UseTpetra){
    const xt_vector_t *tvector = dynamic_cast<const xt_vector_t *>(vec);

    if (tvector->getLocalLength() > 0){
      // getData hangs if vector length is 0
      ArrayRCP<const scalar_t> data = tvector->getData(0);
      elements = data.get();
    }
  }
  else if (map_->lib() == Xpetra::UseEpetra){
    const xe_vector_t *evector = dynamic_cast<const xe_vector_t *>(vec);
      
    if (evector->getLocalLength() > 0){
      // getData hangs if vector length is 0
      ArrayRCP<const double> data = evector->getData(0);

      // Cast so this will compile when scalar_t is not double,
      // a case when this code should never execute.
      elements = reinterpret_cast<const scalar_t *>(data.get());
    }
  }
  else{
    throw std::logic_error("invalid underlying lib");
  }

  ArrayView<const gid_t> gids = map_->getNodeElementList();
  Ids = gids.getRawPtr();

  return getLocalLength();
}

template <typename User>
  template <typename User2>
    size_t XpetraVectorInput<User>::applyPartitioningSolution(
      const User &in, User *&out, 
      const PartitioningSolution<User2> &solution) const
{ 
  // Get an import list

  size_t len = solution.getNumberOfIds();
  const gid_t *gids = solution.getGlobalIdList();
  const size_t *parts = solution.getPartList();
  ArrayRCP<gid_t> gidList = arcp(const_cast<gid_t *>(gids), 0, len, false);
  ArrayRCP<size_t> partList = arcp(const_cast<size_t *>(parts), 0, len, false);

  ArrayRCP<lno_t> dummyIn;
  ArrayRCP<gid_t> importList;
  ArrayRCP<lno_t> dummyOut;
  size_t numNewRows;

  const RCP<const Comm<int> > comm = map_->getComm(); 

  try{
    numNewRows = convertPartListToImportList<gid_t, lno_t, lno_t>(
      *comm, partList, gidList, dummyIn, importList, dummyOut);
  }
  Z2_FORWARD_EXCEPTIONS;

  RCP<const User> inPtr = rcp(&in, false);
  lno_t localNumElts = numNewRows;

  RCP<const User> outPtr = XpetraTraits<User>::doMigration(
   inPtr, localNumElts, importList.get());

  out = const_cast<User *>(outPtr.get());
  outPtr.release();
  return numNewRows;
}

}  //namespace Zoltan2
  
#endif
