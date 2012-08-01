/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK2_LOCALFILTER_DECL_HPP
#define IFPACK2_LOCALFILTER_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ScalarTraits.hpp"


namespace Ifpack2 {
//! Ifpack2_LocalFilter a class for light-weight extraction of the submatrix corresponding to local rows and columns.

/*! Class Ifpack2_LocalFilter enables a light-weight contruction of an
 Tpetra::RowMatrix-derived object, containing only the elements of the original, 
 distributed matrix with local row and column ID. The local
 submatrix is based on a communicator containing the local process only. 
 Each process will have its local object, corresponding to the local submatrix.
 Submatrices may or may not overlap.
 
 The following instructions can be used to create "localized" matrices:
 \code
 #include "Ifpack2_LocalFilter.hpp"
 ...
 Teuchos::RCP<Tpetra::RowMatrix> A;             // fill the elements of A,
 A->FillComplete();

 Ifpack2_LocalFilter LocalA(A);
 \endcode

 Once created, \c LocalA defined, on each process, the submatrix 
 corresponding to local rows and columns only. The creation 
 and use of
 \c LocalA is "cheap", as the elements of the local matrix are
 obtained through calls to ExtractMyRowCopy on the original, distributed
 matrix, say A. This means that \c A must remain in scope every time 
 \c LocalA is accessed.

 A very convenient use of this class is to use Ifpack2 solvers to
 compute the LU factorizations of local blocks. If applied to
 a localized matrix, several Ifpack2 objects can operator in the same
 phase in a safe way, without non-required data exchange.

 \date 1-Aug-12
 
 */ 

template<class MatrixType>
class  LocalFilter : virtual public Tpetra::RowMatrix<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {
  
public:
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  

public:
  //! Constructor.
  explicit LocalFilter(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& Matrix);
  
  //! Destructor.
  virtual ~LocalFilter();

  //! Returns the communicator.
  virtual const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

  //! Returns the underlying node.
  virtual Teuchos::RCP<Node> getNode() const;

  //! Returns the Map that describes the row distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRowMap() const;

  //! \brief Returns the Map that describes the column distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getColMap() const;

  //! Returns the Map that describes the domain distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;

  //! \brief Returns the Map that describes the range distribution in this matrix.
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;

  //! Returns the RowGraph associated with this matrix. 
  virtual Teuchos::RCP<const Tpetra::RowGraph<LocalOrdinal,GlobalOrdinal,Node> > getGraph() const;

  //! Returns the number of global rows in this matrix.
  virtual global_size_t getGlobalNumRows() const;

  //! \brief Returns the number of global columns in this matrix.
  virtual global_size_t getGlobalNumCols() const;

  //! Returns the number of rows owned on the calling node.
  virtual size_t getNodeNumRows() const;

  //! Returns the number of columns needed to apply the forward operator on this node, i.e., the number of elements listed in the column map.
  virtual size_t getNodeNumCols() const;

  //! Returns the index base for global indices for this matrix. 
  virtual GlobalOrdinal getIndexBase() const;

  //! Returns the global number of entries in this matrix.
  virtual global_size_t getGlobalNumEntries() const;

  //! Returns the local number of entries in this matrix.
  virtual size_t getNodeNumEntries() const;
  
  //! \brief Returns the current number of entries on this node in the specified global row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified global row does not belong to this graph. */
  virtual size_t getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const;
  
  //! Returns the current number of entries on this node in the specified local row.
  /*! Returns Teuchos::OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this graph. */
  virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const;
  
  //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
  virtual global_size_t getGlobalNumDiags() const;
  
  //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
  virtual size_t getNodeNumDiags() const;
  
  //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
  virtual size_t getGlobalMaxNumRowEntries() const;
  
  //! \brief Returns the maximum number of entries across all rows/columns on this node.
  virtual size_t getNodeMaxNumRowEntries() const;
  
  //! \brief Indicates whether this matrix has a well-defined column map. 
  virtual bool hasColMap() const;
  
  //! \brief Indicates whether this matrix is lower triangular.
  virtual bool isLowerTriangular() const;
  
  //! \brief Indicates whether this matrix is upper triangular.
  virtual bool isUpperTriangular() const;
  
  //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
  virtual bool isLocallyIndexed() const;
  
  //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
  virtual bool isGloballyIndexed() const;
  
  //! Returns \c true if fillComplete() has been called.
  virtual bool isFillComplete() const;
  

  //@}
  
  //! @name Extraction Methods
  //@{ 
  
  //! Extract a list of entries in a specified global row of this matrix. Put into pre-allocated storage.
  /*!
    \param LocalRow - (In) Global row number for which indices are desired.
    \param Indices - (Out) Global column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumEntries - (Out) Number of indices.
    
    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c GlobalRow. If \c GlobalRow does not belong to this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void getGlobalRowCopy(GlobalOrdinal GlobalRow,
				const Teuchos::ArrayView<GlobalOrdinal> &Indices,
				const Teuchos::ArrayView<Scalar> &Values,
				size_t &NumEntries) const;
  
  //! Extract a list of entries in a specified local row of the graph. Put into storage allocated by calling routine.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices - (Out) Local column indices corresponding to values.
    \param Values - (Out) Matrix values.
    \param NumIndices - (Out) Number of indices.
    
    Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
    with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
    returned as Teuchos::OrdinalTraits<size_t>::invalid().
  */
  virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
			       const Teuchos::ArrayView<LocalOrdinal> &Indices, 
			       const Teuchos::ArrayView<Scalar> &Values,
			       size_t &NumEntries) const ;
  
  //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
  /*!
    \param GlobalRow - (In) Global row number for which indices are desired.
    \param Indices   - (Out) Global column indices corresponding to values.
    \param Values    - (Out) Row values
    \pre <tt>isLocallyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>
    
    Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getGlobalRowView(GlobalOrdinal GlobalRow, 
				Teuchos::ArrayView<const GlobalOrdinal> &indices, 
				Teuchos::ArrayView<const Scalar> &values) const;
  
  //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
  /*!
    \param LocalRow - (In) Local row number for which indices are desired.
    \param Indices  - (Out) Global column indices corresponding to values.
    \param Values   - (Out) Row values
    \pre <tt>isGloballyIndexed() == false</tt>
    \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>
    
    Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
  */
  virtual void getLocalRowView(LocalOrdinal LocalRow, 
			       Teuchos::ArrayView<const LocalOrdinal> &indices, 
			       Teuchos::ArrayView<const Scalar> &values) const;
  
  //! \brief Get a copy of the diagonal entries owned by this node, with local row indices.
  /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
    the zero and non-zero diagonals owned by this node. */
  virtual void getLocalDiagCopy(Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const;
  
  //@}
  
  //! \name Mathematical Methods
  //@{
  
  /**
   * \brief Scales the RowMatrix on the left with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(i)*A(i,j) 
   * where i denoes the global row number of A and 
   * j denotes the global column number of A.
   *
   * \param x A vector to left scale this matrix.
   */
  virtual void leftScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);
  
  /**
   * \brief Scales the RowMatrix on the right with the Vector x.
   *
   * This matrix will be scaled such that A(i,j) = x(j)*A(i,j) 
   * where i denoes the global row number of A and 
   * j denotes the global column number of A.
   *
   * \param x A vector to right scale this matrix.
   */
  virtual void rightScale(const Tpetra::Vector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x);
  
  //! Returns the Frobenius norm of the matrix. 
  /** Computes and returns the Frobenius norm of the matrix, defined as:
      \f$ \|A\|_F = \sqrt{\sum_{i,j} \|\a_{ij}\|^2} \f$
  */
  virtual typename Teuchos::ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const;

  //! \brief Computes the operator-multivector application.
  /*! Loosely, performs \f$Y = \alpha \cdot A^{\textrm{mode}} \cdot X + \beta \cdot Y\f$. However, the details of operation
    vary according to the values of \c alpha and \c beta. Specifically
    - if <tt>beta == 0</tt>, apply() <b>must</b> overwrite \c Y, so that any values in \c Y (including NaNs) are ignored.
    - if <tt>alpha == 0</tt>, apply() <b>may</b> short-circuit the operator, so that any values in \c X (including NaNs) are ignored.
  */
  virtual void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
		     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
		     Teuchos::ETransp mode = Teuchos::NO_TRANS, 
		     Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
		     Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;
  
  //! Indicates whether this operator supports applying the adjoint operator.
  virtual bool hasTransposeApply() const;
  
  //@}
  
  //! \name Deprecated routines to be removed at some point in the future.
  //@{
  
  //! Deprecated. Get a persisting const view of the entries in a specified global row of this matrix.
  /*!
    \param GlobalRow - (In) Global row from which to retrieve matrix entries.
    \param Indices - (Out) Indices for the global row.
    \param Values - (Out) Values for the global row.
    
    Note: If \c GlobalRow does not belong to this node, then \c Indices and \c Values are set to <tt>Teuchos::null</t>>.
    
    \pre isLocallyIndexed()==false
  */
  TPETRA_DEPRECATED virtual void getGlobalRowView(GlobalOrdinal GlobalRow, 
						  Teuchos::ArrayRCP<const GlobalOrdinal> &indices,
						  Teuchos::ArrayRCP<const Scalar>        &values) const;
  
  //! Deprecated. Get a persisting const view of the entries in a specified local row of this matrix.
  /*!
    \param LocalRow - (In) Local row from which to retrieve matrix entries.
    \param Indices - (Out) Indices for the local row.
    \param Values - (Out) Values for the local row.
    
    Note: If \c LocalRow is not valid for this node, then \c Indices and \c Values are set to <tt>Teuchos::null</tt>.
    
    \pre isGloballyIndexed()==false
  */
  TPETRA_DEPRECATED virtual void getLocalRowView(LocalOrdinal LocalRow,
						 Teuchos::ArrayRCP<const LocalOrdinal> &indices,
						 Teuchos::ArrayRCP<const Scalar>       &values) const;  
  //@}

  
private:

  //! Pointer to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > A_;
  //! Communicator containing this process only.
  Teuchos::RCP<Teuchos::Comm<int> > SerialComm_;
  //! Map based on SerialComm_, containing the local rows only.
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > LocalMap_;
  //! Number of rows in the local matrix.
  size_t NumRows_;
  //! Number of nonzeros in the local matrix.
  size_t NumNonzeros_;
  //! Maximum number of nonzero entries in a row for the filtered matrix.
  size_t MaxNumEntries_;
  //! Maximum number of nonzero entries in a row for the unfiltered matrix.
  size_t MaxNumEntriesA_;
  //! NumEntries_[i] contains the nonzero entries in row `i'.
  std::vector<size_t> NumEntries_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable Teuchos::Array<LocalOrdinal> Indices_;
  //! Used in ExtractMyRowCopy, to avoid allocation each time.
  mutable Teuchos::Array<Scalar> Values_;

  //! Norms
  Scalar NormOne_;
  Scalar NormInf_;

};// class LocalFilter

}// namespace Ifpack2

#endif /* IFPACK2_LOCALFILTER_DECL_HPP */
