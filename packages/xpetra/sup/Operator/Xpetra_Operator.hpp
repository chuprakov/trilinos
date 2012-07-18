#ifndef XPETRA_OPERATOR_HPP
#define XPETRA_OPERATOR_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_Exceptions.hpp"

#include "Xpetra_MultiVector.hpp"
#include "Xpetra_CrsGraph.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_CrsMatrixFactory.hpp"
#include "Xpetra_OperatorView.hpp"
#include "Xpetra_StridedMap.hpp"
#include "Xpetra_StridedMapFactory.hpp"

#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_Hashtable.hpp>

/** \file Xpetra_Operator.hpp

Declarations for the class Xpetra::Operator.
*/
namespace Xpetra {

  typedef std::string viewLabel_t;

  template <class Scalar, 
            class LocalOrdinal  = int, 
            class GlobalOrdinal = LocalOrdinal, 
            class Node          = Kokkos::DefaultNode::DefaultNodeType, 
            class LocalMatOps   = typename Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Operator : virtual public Teuchos::Describable {
    
    typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrix;
    typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#ifdef HAVE_XPETRA_TPETRA
    typedef Xpetra::TpetraCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsMatrix;
#endif
    typedef Xpetra::CrsMatrixFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixFactory;
    typedef Xpetra::OperatorView<LocalOrdinal, GlobalOrdinal, Node> OperatorView;

  public:
  
    //! @name Constructor/Destructor Methods
    //@{

    Operator() : blksize_(1) { }

    //! Destructor
    virtual ~Operator() { }

    //@}

    //! @name View management methods
    //@{
    void CreateView(viewLabel_t viewLabel, const RCP<const Map> & rowMap, const RCP<const Map> & colMap) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == true, Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.CreateView(): a view labeled '" + viewLabel + "' already exist.");
      RCP<OperatorView> view = rcp(new OperatorView(rowMap, colMap));
      operatorViewTable_.put(viewLabel, view);
    }
    
    void CreateView(const viewLabel_t viewLabel, const RCP<Operator> & A, bool transposeA = false, const RCP<Operator> & B = Teuchos::null, bool transposeB = false) {
            
      RCP<const Map> domainMap = Teuchos::null;
      RCP<const Map> rangeMap  = Teuchos::null;

      // A has strided Maps
      if(A->IsView(viewLabel)) {
        Xpetra::viewLabel_t oldView = A->SwitchToView(viewLabel); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
        rangeMap = (transposeA) ? A->getColMap() : A->getRowMap();
        domainMap = (transposeA) ? A->getRowMap() : A->getColMap(); // overwrite if B != Teuchos::null
        oldView = A->SwitchToView(oldView);
      } else rangeMap = (transposeA) ? A->getDomainMap() : A->getRangeMap();
      
      // B has strided Maps
      if(B != Teuchos::null ) {
        if(B->IsView(viewLabel)) {
          Xpetra::viewLabel_t oldView = B->SwitchToView(viewLabel); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
          domainMap = (transposeB) ? B->getRowMap() : B->getColMap();
          oldView = B->SwitchToView(oldView);
        } else domainMap = (transposeB) ? B->getRangeMap() : B->getDomainMap();
      }
      
      
      if(IsView(viewLabel)) RemoveView(viewLabel);
      
      CreateView(viewLabel, rangeMap, domainMap);   
    }

    
    void RemoveView(const viewLabel_t viewLabel) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.RemoveView(): view '" + viewLabel + "' does not exist.");
      TEUCHOS_TEST_FOR_EXCEPTION(viewLabel == GetDefaultViewLabel(), Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.RemoveView(): view '" + viewLabel + "' is the default view and cannot be removed.");
      operatorViewTable_.remove(viewLabel);
    }
    
    const viewLabel_t SwitchToView(const viewLabel_t viewLabel) {
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.SwitchToView(): view '" + viewLabel + "' does not exist.");
      viewLabel_t oldViewLabel = GetCurrentViewLabel();
      currentViewLabel_ = viewLabel;
      return oldViewLabel;
    }

    bool IsView(const viewLabel_t viewLabel) const {
      return operatorViewTable_.containsKey(viewLabel);
    }

    const viewLabel_t SwitchToDefaultView() { return SwitchToView(GetDefaultViewLabel()); }

    const viewLabel_t & GetDefaultViewLabel() const { return defaultViewLabel_; }

    const viewLabel_t & GetCurrentViewLabel() const { return currentViewLabel_; }
    
    //@}

    //! @name Insertion/Removal Methods
    //@{ 

    //! Insert matrix entries, using global IDs.
    /** All index values must be in the global space. 
        \pre \c globalRow exists as an ID in the global row map
        \pre <tt>isLocallyIndexed() == false</tt>
        \pre <tt>isStorageOptimized() == false</tt>

        \post <tt>isGloballyIndexed() == true</tt>

        \note If \c globalRow does not belong to the matrix on this node, then it will be communicated to the appropriate node when globalAssemble() is called (which will, at the latest, occur during the next call to fillComplete().) Otherwise, the entries will be inserted in the local matrix. 
        \note If the matrix row already contains values at the indices corresponding to values in \c cols, then the new values will be summed with the old values; this may happen at insertion or during the next call to fillComplete().
        \note If <tt>hasColMap() == true</tt>, only (cols[i],vals[i]) where cols[i] belongs to the column map on this node will be inserted into the matrix.
    */
    virtual void insertGlobalValues(GlobalOrdinal globalRow, const ArrayView<const GlobalOrdinal> &cols, const ArrayView<const Scalar> &vals) = 0;

    //@}

    //! @name Transformational Methods
    //@{ 

    /*! \brief Signal that data entry is complete, specifying domain and range maps.

    Off-node indices are distributed (via globalAssemble()), indices are sorted, redundant indices are eliminated, and global indices are transformed to local indices.

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */ 
    virtual void fillComplete(const RCP<const Map> &domainMap, const RCP<const Map> &rangeMap, const RCP<ParameterList> &params = null) =0;

    /*! \brief Signal that data entry is complete. 

    Off-node entries are distributed (via globalAssemble()), repeated entries are summed, and global indices are transformed to local indices.

    \note This method calls fillComplete( getRowMap(), getRowMap(), os ).

    \pre  <tt>isFillActive() == true<tt>
    \pre <tt>isFillComplete()() == false<tt>

    \post <tt>isFillActive() == false<tt>
    \post <tt>isFillComplete() == true<tt>
    \post if <tt>os == DoOptimizeStorage<tt>, then <tt>isStorageOptimized() == true</tt>
    */
    //TODO : Get ride of "Tpetra"::OptimizeOption
    virtual void fillComplete(const RCP<ParameterList> &params=null) =0;

    //@}

    //! @name Methods implementing RowMatrix
    //@{ 

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP<const Map> & getRowMap() const { return getRowMap(GetCurrentViewLabel()); }

    //! Returns the Map that describes the row distribution in this matrix.
    virtual const RCP<const Map> & getRowMap(viewLabel_t viewLabel) const { 
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.GetRowMap(): view '" + viewLabel + "' does not exist.");
      return operatorViewTable_.get(viewLabel)->GetRowMap(); 
    }

    //! \brief Returns the Map that describes the column distribution in this matrix.
    //! This might be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Map> & getColMap() const { return getColMap(GetCurrentViewLabel()); }

    //! \brief Returns the Map that describes the column distribution in this matrix.
    virtual const RCP<const Map> & getColMap(viewLabel_t viewLabel) const { 
      TEUCHOS_TEST_FOR_EXCEPTION(operatorViewTable_.containsKey(viewLabel) == false, Xpetra::Exceptions::RuntimeError, "Xpetra::Operator.GetColMap(): view '" + viewLabel + "' does not exist.");
      return operatorViewTable_.get(viewLabel)->GetColMap(); 
    }

    //! Returns the number of global rows in this matrix.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumRows() const =0;

    //! \brief Returns the number of global columns in the matrix.
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumCols() const =0;

    //! Returns the number of matrix rows owned on the calling node.
    virtual size_t getNodeNumRows() const =0;

    //! Returns the global number of entries in this matrix.
    virtual global_size_t getGlobalNumEntries() const =0;

    //! Returns the local number of entries in this matrix.
    virtual size_t getNodeNumEntries() const =0;

    //! Returns the current number of entries on this node in the specified local row.
    /*! Returns OrdinalTraits<size_t>::invalid() if the specified local row is not valid for this matrix. */
    virtual size_t getNumEntriesInLocalRow(LocalOrdinal localRow) const =0;

    //! \brief Returns the number of global diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual global_size_t getGlobalNumDiags() const =0;

    //! \brief Returns the number of local diagonal entries, based on global row/column index comparisons. 
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeNumDiags() const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on all nodes.
    /** Undefined if isFillActive().
     */
    virtual size_t getGlobalMaxNumRowEntries() const =0;

    //! \brief Returns the maximum number of entries across all rows/columns on this node.
    /** Undefined if isFillActive().
     */
    virtual size_t getNodeMaxNumRowEntries() const =0;

    //! \brief If matrix indices are in the local range, this function returns true. Otherwise, this function returns false. */
    virtual bool isLocallyIndexed() const =0;

    //! \brief If matrix indices are in the global range, this function returns true. Otherwise, this function returns false. */
    virtual bool isGloballyIndexed() const =0;

    //! Returns \c true if fillComplete() has been called and the matrix is in compute mode.
    virtual bool isFillComplete() const =0;

    //! Extract a list of entries in a specified local row of the matrix. Put into storage allocated by calling routine.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices - (Out) Local column indices corresponding to values.
      \param Values - (Out) Matrix values.
      \param NumIndices - (Out) Number of indices.

      Note: A std::runtime_error exception is thrown if either \c Indices or \c Values is not large enough to hold the data associated
      with row \c LocalRow. If \c LocalRow is not valid for this node, then \c Indices and \c Values are unchanged and \c NumIndices is 
      returned as OrdinalTraits<size_t>::invalid().

      \pre <tt>isLocallyIndexed()==true</tt> or <tt>hasColMap() == true</tt>
    */
    virtual void getLocalRowCopy(LocalOrdinal LocalRow, 
                                 const ArrayView<LocalOrdinal> &Indices, 
                                 const ArrayView<Scalar> &Values,
                                 size_t &NumEntries
                                 ) const =0;

    //! Extract a const, non-persisting view of global indices in a specified row of the matrix.
    /*!
      \param GlobalRow - (In) Global row number for which indices are desired.
      \param Indices   - (Out) Global column indices corresponding to values.
      \param Values    - (Out) Row values
      \pre <tt>isLocallyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInGlobalRow(GlobalRow)</tt>

      Note: If \c GlobalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getGlobalRowView(GlobalOrdinal GlobalRow, ArrayView<const GlobalOrdinal> &indices, ArrayView<const Scalar> &values) const =0;

    //! Extract a const, non-persisting view of local indices in a specified row of the matrix.
    /*!
      \param LocalRow - (In) Local row number for which indices are desired.
      \param Indices  - (Out) Global column indices corresponding to values.
      \param Values   - (Out) Row values
      \pre <tt>isGloballyIndexed() == false</tt>
      \post <tt>indices.size() == getNumEntriesInLocalRow(LocalRow)</tt>

      Note: If \c LocalRow does not belong to this node, then \c indices is set to null.
    */
    virtual void getLocalRowView(LocalOrdinal LocalRow, ArrayView<const LocalOrdinal> &indices, ArrayView<const Scalar> &values) const =0;

    //! \brief Get a copy of the diagonal entries owned by this node, with local row idices.
    /*! Returns a distributed Vector object partitioned according to this matrix's row map, containing the 
      the zero and non-zero diagonals owned by this node. */
    virtual void getLocalDiagCopy(Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &diag) const =0;

    //! Get Frobenius norm of the matrix
    virtual typename ScalarTraits<Scalar>::magnitudeType getFrobeniusNorm() const = 0;

    //@}

    //! @name Advanced Matrix-vector multiplication and solve methods
    //@{

    //! Multiplies this matrix by a MultiVector.
    /*! \c X is required to be post-imported, i.e., described by the column map of the matrix. \c Y is required to be pre-exported, i.e., described by the row map of the matrix.

    Both are required to have constant stride, and they are not permitted to ocupy overlapping space. No runtime checking will be performed in a non-debug build.

    This method is templated on the scalar type of MultiVector objects, allowing this method to be applied to MultiVector objects of arbitrary type. However, it is recommended that multiply() not be called directly; instead, use the CrsMatrixMultiplyOp, as it will handle the import/exprt operations required to apply a matrix with non-trivial communication needs.
          
    If \c beta is equal to zero, the operation will enjoy overwrite semantics (\c Y will be overwritten with the result of the multiplication). Otherwise, the result of the multiplication
    will be accumulated into \c Y.
    */
    //TODO virtual=0 // TODO: Add default parameters ?
//     virtual void multiply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, Teuchos::ETransp trans, Scalar alpha, Scalar beta) const=0;

    //@}

    //! @name Methods implementing Operator
    //@{ 

    //! \brief Computes the sparse matrix-multivector multiplication.
    /*! Performs \f$Y = \alpha A^{\textrm{mode}} X + \beta Y\f$, with one special exceptions:
      - if <tt>beta == 0</tt>, apply() overwrites \c Y, so that any values in \c Y (including NaNs) are ignored.
    */
    virtual void apply(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & X, MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
                       Teuchos::ETransp mode = Teuchos::NO_TRANS,
                       Scalar alpha = ScalarTraits<Scalar>::one(),
                       Scalar beta = ScalarTraits<Scalar>::zero()) const =0;

    //! \brief Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getDomainMap() const =0;

    //! Returns the Map associated with the domain of this operator.
    //! This will be <tt>null</tt> until fillComplete() is called.
    virtual const RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > getRangeMap() const =0;

    //@}

    //! @name Overridden from Teuchos::Describable 
    //@{

    // TODO: describe of views can be done here

    //   /** \brief Return a simple one-line description of this object. */
    //   virtual std::string description() const =0;

    //   /** \brief Print the object with some verbosity level to an FancyOStream object. */
    //   virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;

    //@}


    //! @name Overridden from Teuchos::Describable
    //@{

    /** \brief Return a simple one-line description of this object. */
    virtual std::string description() const =0;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;
    //@}

    // JG: Added:

    //! Returns the CrsGraph associated with this matrix. 
    virtual RCP<const CrsGraph> getCrsGraph() const =0;

  public: //TODO: protected
    Teuchos::Hashtable<viewLabel_t, RCP<OperatorView> > operatorViewTable_; // hashtable storing the operator views (keys = view names, values = views).
  
    viewLabel_t defaultViewLabel_;  // label of the view associated with inital Operator construction
    viewLabel_t currentViewLabel_;  // label of the current view

    // ----------------------------------------------------------------------------------
    // "TEMPORARY" VIEW MECHANISM
    // TODO: the view mechanism should be implemented as in MueMat.
    LocalOrdinal blksize_; // TODO remove this
    // RCP<GOVector> variableBlockSizeInfo_; TODO: should be moved from CoalesceDropFactory to here.

    void SetFixedBlockSize(LocalOrdinal blksize) { 
      TEUCHOS_TEST_FOR_EXCEPTION(isFillComplete() == false, Exceptions::RuntimeError, "Xpetra::Operator::SetFixedBlockSize(): operator is not filled and completed."); // TODO: do we need this? we just wanna "copy" the domain and range map

      std::vector<size_t> stridingInfo;
      stridingInfo.push_back(Teuchos::as<size_t>(blksize));
      LocalOrdinal stridedBlockId = -1;
            
      RCP<const Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> > stridedRangeMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(getRangeMap()->lib(),
                                                    getRangeMap(),
                                                    stridingInfo,
                                                    stridedBlockId
                                                    );
      RCP<const Map> stridedDomainMap = Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node>::Build(getDomainMap()->lib(),
					      getDomainMap(),
					      stridingInfo,
					      stridedBlockId
					      );

      if(IsView("stridedMaps") == true) RemoveView("stridedMaps");
      CreateView("stridedMaps", stridedRangeMap, stridedDomainMap); 
    }
    LocalOrdinal GetFixedBlockSize() const {
      if(IsView("StridedMaps")==true) {
        Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> > rangeMap = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >(getRowMap("stridedMaps"));
        Teuchos::RCP<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> > domainMap = Teuchos::rcp_dynamic_cast<const StridedMap<LocalOrdinal, GlobalOrdinal, Node> >(getColMap("stridedMaps"));
        TEUCHOS_TEST_FOR_EXCEPTION(rangeMap  == Teuchos::null, Exceptions::BadCast, "Xpetra::Operator::GetFixedBlockSize(): rangeMap is not of type StridedMap");
        TEUCHOS_TEST_FOR_EXCEPTION(domainMap == Teuchos::null, Exceptions::BadCast, "Xpetra::Operator::GetFixedBlockSize(): domainMap is not of type StridedMap");
        TEUCHOS_TEST_FOR_EXCEPTION(domainMap->getFixedBlockSize() != rangeMap->getFixedBlockSize(), Exceptions::RuntimeError, "Xpetra::Operator::GetFixedBlockSize(): block size of rangeMap and domainMap are different.");
        TEUCHOS_TEST_FOR_EXCEPTION(domainMap->getFixedBlockSize() != Teuchos::as<size_t>(blksize_), Exceptions::RuntimeError, "Xpetra::Operator::GetFixedBlockSize(): block size of domainMap and blksize_ are different.");
        return Teuchos::as<LocalOrdinal>(domainMap->getFixedBlockSize()); // TODO: why LocalOrdinal?
      } else
        TEUCHOS_TEST_FOR_EXCEPTION(false, Exceptions::RuntimeError, "Xpetra::Operator::GetFixedBlockSize(): no blksize_ available."); // TODO remove this
        return blksize_;
    }; //TODO: why LocalOrdinal?

    // ----------------------------------------------------------------------------------

  }; //class Operator

} //namespace Xpetra

#define XPETRA_OPERATOR_SHORT
#endif //XPETRA_OPERATOR_DECL_HPP
