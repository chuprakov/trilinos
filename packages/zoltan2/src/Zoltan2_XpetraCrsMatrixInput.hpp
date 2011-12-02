// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_XpetraCrsMatrixInput.hpp

    \brief An input adapter for a Xpetra::CrsMatrix.
*/

#ifndef _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_
#define _ZOLTAN2_XPETRACRSMATRIXINPUT_HPP_

#include <Xpetra_CrsMatrix.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Zoltan2_MatrixInput.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_Util.hpp>

namespace Zoltan2 {

//////////////////////////////////////////////////////////////////////////////
/*! Zoltan2::XpetraCrsMatrixInput
    \brief Provides access for Zoltan2 to Xpetra::CrsMatrix data.

    TODO: we assume FillComplete has been called.  We should support
                objects that are not FillCompleted.

    The template parameter is the user's input object - an Epetra
    matrix or a templated Tpetra::CrsMatrix 
    or a templated Xpetra::CrsMatrix.
    TODO add RowMatrix
*/

template <typename User>
class XpetraCrsMatrixInput : public MatrixInput<User> {
public:

  typedef typename InputTraits<User>::scalar_t scalar_t;
  typedef typename InputTraits<User>::lno_t    lno_t;
  typedef typename InputTraits<User>::gno_t    gno_t;
  typedef typename InputTraits<User>::lid_t    lid_t;
  typedef typename InputTraits<User>::gid_t    gid_t;
  typedef typename InputTraits<User>::node_t   node_t;
  typedef Xpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> xmatrix_t;

  /*! Destructor
   */
  ~XpetraCrsMatrixInput() { }

  /*! Constructor   
   */
  // Constructor 
  XpetraCrsMatrixInput(const RCP<const User> &inmatrix):
    inmatrix_(inmatrix), 
    matrix_(),
    rowMap_(),
    colMap_(),
    base_(),
    offset_(),
    columnIds_()
  {
    matrix_ = XpetraTraits<User>::convertToXpetra(inmatrix);
    rowMap_ = matrix_->getRowMap();
    colMap_ = matrix_->getColMap();
    base_ = rowMap_->getIndexBase();

    size_t nrows = matrix_->getNodeNumRows();
    size_t nnz = matrix_->getNodeNumEntries();
 
    offset_.resize(nrows+1, lid_t(0));
    columnIds_.resize(nnz);
    ArrayView<const lid_t> indices;
    ArrayView<const scalar_t> nzs;
    lid_t next = 0;
    for (size_t i=0; i < nrows; i++){
      lid_t row = i + base_;
      lid_t nnz = matrix_->getNumEntriesInLocalRow(row);
      matrix_->getLocalRowView(row, indices, nzs);
      for (lid_t j=0; j < nnz; j++){
        // TODO - this will be slow
        //   Is it possible that global columns ids might be stored in order?
        columnIds_[next++] = colMap_->getGlobalElement(indices[j]);
      }
      offset_[i+1] = offset_[i] + nnz;
    } 
  };

  /*! Access to xpetra matrix
   */

  const RCP<const xmatrix_t> &getMatrix() const
  {
    return matrix_;
  }

  ////////////////////////////////////////////////////
  // The InputAdapter interface.
  ////////////////////////////////////////////////////

  std::string inputAdapterName()const {
    return std::string("XpetraCrsMatrix");}

  bool haveLocalIds() const { return true;}

  bool haveConsecutiveLocalIds(size_t &base) const{
    base = base_;
    return true;
  }

  ////////////////////////////////////////////////////
  // The MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Returns the number rows on this process.
   */
  size_t getLocalNumRows() const { 
    return matrix_->getNodeNumRows();
  }

  /*! Returns the number rows in the entire matrix.
   */
  global_size_t getGlobalNumRows() const { 
    return matrix_->getGlobalNumRows();
  }

  /*! Returns the number columns on this process.
   */
  size_t getLocalNumColumns() const { 
    return matrix_->getNodeNumCols();
  }

  /*! Returns the number columns on this entire matrix.
   *    what about directional columns, count twice?
   */
  global_size_t getGlobalNumColumns() const { 
    return matrix_->getGlobalNumCols();
  }

  /*! Return a read only view of the data.
     \param rowIds  Global row ids.  The memory for the global 
          row IDs persists until the underlying Xpetra::CrsMatrix is deleted.
     \param localIds on return is NULL, signifying that local IDs are
           contiguous integers starting at the base supplied in
          haveConsecutiveLocalIds.
     \param offsets The columns for rowIds[i] begin at colIds[offsets[i]].  
        There are numRows+1 offsets.  The last offset is the length of the 
        colIds array.  The memory pointed to by offsets persists 
        in the GraphModel is deleted.
     \param colIds The global column Ids. The memory pointed to by colIds 
          persists until the GraphModel is deleted.

     \return  The number rows in the rowIds list is returned.
   */

  size_t getRowListView(const gid_t *&rowIds, const lid_t *&localIds,
    const lid_t *&offsets, const gid_t *& colIds) const
  {
    size_t nrows = getLocalNumRows();

    ArrayView<const gid_t> rowView = rowMap_->getNodeElementList();
    rowIds = rowView.getRawPtr();
   
    localIds = NULL;   // Implies consecutive integers

    offsets = offset_.getRawPtr();
    colIds = columnIds_.getRawPtr();
    return nrows;
  }

  ////////////////////////////////////////////////////
  // End of MatrixInput interface.
  ////////////////////////////////////////////////////

  /*! Apply a partitioning solution to the matrix.
   *   Every gid that was belongs to this process must
   *   be on the list, or the Import will fail.
   *
   *   TODO : params etc
   */
  size_t applyPartitioningSolution(const User &in, User *&out,
         const PartitioningSolution<gid_t, lid_t, lno_t> &solution)
  { 
    // Get an import list

    ArrayRCP<gid_t>  gidList  = solution.getGidsRCPConst();
    ArrayRCP<size_t> partList = solution.getPartsRCPConst();
    ArrayRCP<lno_t> dummyIn;
    ArrayRCP<gid_t> importList;
    ArrayRCP<lno_t> dummyOut;
    size_t numNewRows;
    const RCP<const Comm<int> > comm = matrix_->getRowMap()->getComm();

    try{
      numNewRows = convertPartListToImportList<gid_t, lno_t, lno_t>(
        *comm, partList, gidList, dummyIn, importList, dummyOut);
    }
    Z2_FORWARD_EXCEPTIONS;

    gno_t lsum = numNewRows;
    gno_t gsum = 0;
    Teuchos::reduceAll<int, gno_t>(*comm, Teuchos::REDUCE_SUM, 1, &lsum, &gsum);

    RCP<const User> inPtr = rcp(&in, false);

    RCP<const User> outPtr = XpetraTraits<User>::doMigration(
     inPtr, lsum, importList.getRawPtr(), base_);

    out = const_cast<User *>(outPtr.get());
    outPtr.release();
    return numNewRows;
  }

private:

  RCP<const User> inmatrix_;
  RCP<const xmatrix_t> matrix_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > rowMap_;
  RCP<const Xpetra::Map<lno_t, gno_t, node_t> > colMap_;
  lno_t base_;
  ArrayRCP<lno_t> offset_;
  ArrayRCP<gno_t> columnIds_;
};
  
}  //namespace Zoltan2
  
#endif
