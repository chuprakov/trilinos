#ifndef CTHULHU_TPETRABLOCKMAP_HPP
#define CTHULHU_TPETRABLOCKMAP_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_TPETRA
#error This file should be included only if HAVE_CTHULHU_TPETRA is defined.
#endif

#include "Cthulhu_Map.hpp"
#include "Cthulhu_BlockMap.hpp"

#include "Tpetra_BlockMap.hpp"
#include "Tpetra_Map.hpp"

#include "Cthulhu_Debug.hpp"
#include "Cthulhu_Exceptions.hpp"

/** \file Cthulhu_EpetraBlockMap.hpp

  Declarations for the class Cthulhu::TpetraBlockMap.
*/
namespace Cthulhu {

/** \brief Block-entry counterpart to Cthulhu::Map.

  BlockMap doesn't inherit Cthulhu::Map
*/
template <class LocalOrdinal, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType>
class TpetraBlockMap : public Cthulhu::BlockMap<LocalOrdinal,GlobalOrdinal,Node> {
 public:
  
  // The following typedef are used by the CTHULHU_DYNAMIC_CAST() macro.
  typedef TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;

  //! @name Constructor/Destructor Methods
  //@{

  /*! \brief TpetraBlockMap constructor specifying numGlobalBlocks and constant blockSize.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 LocalOrdinal blockSize,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, blockSize, indexBase, comm, node))) { CTHULHU_DEBUG_ME; }

  /*! \brief TpetraBlockMap constructor specifying num global and local blocks, and constant blockSize.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 size_t numLocalBlocks,
                 LocalOrdinal blockSize,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, numLocalBlocks, blockSize, indexBase, comm, node))) { CTHULHU_DEBUG_ME; }

  /*! \brief TpetraBlockMap constructor specifying numGlobalBlocks and lists of local blocks first-global-point-in-blocks, and blockSizes.
   */
  TpetraBlockMap(global_size_t numGlobalBlocks,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myFirstGlobalPointInBlocks,
                 const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
                 GlobalOrdinal indexBase,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) : map_(rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(numGlobalBlocks, myGlobalBlockIDs, myFirstGlobalPointInBlocks, myBlockSizes, indexBase, comm, node))) { CTHULHU_DEBUG_ME; }
  
  /*! \brief TpetraBlockMap constructor which takes a "regular" Map.
   * The arrays myGlobalBlockIDs and myBlockSizes must be the same length, and
   * sum(myBlockSizes) must equal pointMap->getNodeNumElements().
   * If these arrays are different lengths or sum(myBlockSizes) is incorrect,
   * then std::runtime_error is thrown.
   */
  TpetraBlockMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& pointMap,
                 const Teuchos::ArrayView<const GlobalOrdinal>& myGlobalBlockIDs,
                 const Teuchos::ArrayView<const LocalOrdinal>& myBlockSizes,
                 const Teuchos::RCP<Node> &node = Kokkos::DefaultNode::getDefaultNode()) {
    CTHULHU_DEBUG_ME;
    CTHULHU_RCP_DYNAMIC_CAST(const TpetraMap, pointMap, tPointMap, "Cthulhu::TpetraBlockMap constructors only accept Cthulhu::TpetraMap as input arguments.");
    map_ = rcp(new Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node>(tPointMap->getTpetra_Map(), myGlobalBlockIDs, myBlockSizes, node));
  }

  TpetraBlockMap(const Teuchos::RCP<const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > &map) : map_(map) { CTHULHU_DEBUG_ME; }

  //! TpetraBlockMap destructor.
  virtual ~TpetraBlockMap(){ CTHULHU_DEBUG_ME; }

  //@}

  //! @name Attribute Accessor Methods
  //@{
#ifdef CTHULHU_NOT_IMPLEMENTED
  inline const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> >& getPointMap() const { CTHULHU_DEBUG_ME; return map_->getPointMap(); }
#endif

  inline global_size_t getGlobalNumBlocks() const { CTHULHU_DEBUG_ME; return map_->getGlobalNumBlocks(); }

  //! Return number of blocks on the local processor.
  inline size_t getNodeNumBlocks() const { CTHULHU_DEBUG_ME; return map_->getNodeNumBlocks(); }

  inline Teuchos::ArrayView<const GlobalOrdinal> getNodeBlockIDs() const { CTHULHU_DEBUG_ME; return map_->getNodeBlockIDs(); }

  inline bool isBlockSizeConstant() const { CTHULHU_DEBUG_ME; return map_->isBlockSizeConstant(); }

  //! Return ArrayRCP of first-local-point in local blocks.
  inline Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks() const { CTHULHU_DEBUG_ME; return map_->getNodeFirstPointInBlocks(); }

  //! Return device-resident ArrayRCP of first-local-point in local blocks.
  /*! This version of this method is primarily used internally by VbrMatrix
      for passing data to the matrix-vector-product kernel.
  */
  inline Teuchos::ArrayRCP<const LocalOrdinal> getNodeFirstPointInBlocks_Device() const { CTHULHU_DEBUG_ME; return map_->getNodeFirstPointInBlocks_Device(); }

  //! Return the globalBlockID corresponding to the given localBlockID
  /*! If localBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  inline GlobalOrdinal getGlobalBlockID(LocalOrdinal localBlockID) const { CTHULHU_DEBUG_ME; return map_->getGlobalBlockID(localBlockID); }

  //! Return the localBlockID corresponding to the given globalBlockID
  /*! If globalBlockID is not present on this processor, returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid().
  */
  inline LocalOrdinal getLocalBlockID(GlobalOrdinal globalBlockID) const { CTHULHU_DEBUG_ME; return map_->getLocalBlockID(globalBlockID); }

  //! Return the block-size for localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  inline LocalOrdinal getLocalBlockSize(LocalOrdinal localBlockID) const { CTHULHU_DEBUG_ME; return map_->getLocalBlockSize(localBlockID); }

  //! Return the first local point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  inline LocalOrdinal getFirstLocalPointInLocalBlock(LocalOrdinal localBlockID) const { CTHULHU_DEBUG_ME; return map_->getFirstLocalPointInLocalBlock(localBlockID); }

  //! Return the first global point-index corresponding to localBlockID
  /*! If localBlockID is out of range (less than 0 or greater/equal num-local-blocks),
   * then std::runtime_error is thrown.
   */
  inline GlobalOrdinal getFirstGlobalPointInLocalBlock(LocalOrdinal localBlockID) const { CTHULHU_DEBUG_ME; return map_->getFirstGlobalPointInLocalBlock(localBlockID); }

  //@}

  RCP< const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_BlockMap() const { CTHULHU_DEBUG_ME; return map_; }

private:
  RCP< const Tpetra::BlockMap<LocalOrdinal, GlobalOrdinal, Node> > map_;


};//class TpetraBlockMap

// //-----------------------------------------------------------------
// template<class LocalOrdinal,class GlobalOrdinal,class Node>
// Teuchos::RCP<const Cthulhu::Map<LocalOrdinal,GlobalOrdinal,Node> >
// convertTpetraBlockMapToTpetraPointMap(const Teuchos::RCP<const Cthulhu::TpetraBlockMap<LocalOrdinal,GlobalOrdinal,Node> >& blockMap) { CTHULHU_DEBUG_ME;
//   return rcp(new TpetraMap(convertTpetraBlockMapToTpetraPointMap(blockMap.getTpetra_BlockMap())));
// }

}//namespace Cthulhu

#define CTHULHU_TPETRABLOCKMAP_SHORT
#endif
