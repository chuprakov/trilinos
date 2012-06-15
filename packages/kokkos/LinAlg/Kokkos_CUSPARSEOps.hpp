//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef KOKKOS_CUSPARSEOPS_HPP
#define KOKKOS_CUSPARSEOPS_HPP

#include <Teuchos_DataAccess.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include <cusparse_v2.h>

namespace Kokkos {

  namespace CUSPARSEdetails {
    class CUSPARSEDestroyer {
    public:
      CUSPARSEDestroyer();
      void free(void *ptr);
    };
    static RCP<cusparseHandle_t> session_handle;
    void initCUSPARSEsession();
  }

  //! \class CUSPARSECrsGraph
  /** \brief CRS sparse graph class supporting the CUSPARSE library.
  */
  template <class Node>
  class CUSPARSECrsGraph : public CrsGraphBase<int,Node> 
  {
    public:
      CUSPARSECrsGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params);
      bool isEmpty() const;
      void setStructure(const ArrayRCP<const size_t>  &ptrs,
                        const ArrayRCP<const int> &inds);
      inline ArrayRCP<const size_t> getPointers() const;
      inline ArrayRCP<const int>    getIndices() const;
      inline bool isInitialized() const;
    private:
      ArrayRCP<const size_t>  ptrs_;
      ArrayRCP<const int> inds_;
      bool isInitialized_;
      bool isEmpty_;
  };

  //! \class CUSPARSECrsMatrix 
  /** \brief CRS sparse matrix class supporting the CUSPARSE library.
  */
  template <class Scalar, 
            class Node> 
  class CUSPARSECrsMatrix : public CrsMatrixBase<Scalar,int,Node> 
  {
    public:
      CUSPARSECrsMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<ParameterList> &params);
      void setValues(const ArrayRCP<const Scalar> &vals);
      inline ArrayRCP<const Scalar> getValues() const;
      inline bool isInitialized() const;
    private:
      ArrayRCP<const Scalar> vals_;
      bool isInitialized_;
  };

  template <class Node>
  CUSPARSECrsGraph<Node>::CUSPARSECrsGraph(size_t numRows, const RCP<Node> &node, const RCP<ParameterList> &params)
  : CrsGraphBase<int,Node>(numRows,node,params)
  , isInitialized_(false)
  , isEmpty_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isEmpty() const
  { return isEmpty_; }

  template <class Node>
  void CUSPARSECrsGraph<Node>::setStructure(
                      const ArrayRCP<const size_t>  &ptrs,
                      const ArrayRCP<const int> &inds)
  {
    std::string tfecfFuncName("setStructure(ptrs,inds)");
    const size_t numrows = this->getNumRows();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        (size_t)ptrs.size() != numrows+1 
        || ptrs[0] != 0
        || (size_t)inds.size() != ptrs[numrows],
        std::runtime_error, " graph data not coherent."
    )
    const size_t numEntries = ptrs[numrows];
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix has already been initialized"
    )
    if (numrows == 0 || numEntries == 0) isEmpty_ = true;
    ptrs_ = ptrs;
    inds_ = inds;
    isInitialized_ = true;
  }

  template <class Node>
  ArrayRCP<const size_t> CUSPARSECrsGraph<Node>::getPointers() const
  { return ptrs_; }

  template <class Node>
  ArrayRCP<const int> CUSPARSECrsGraph<Node>::getIndices() const
  { return inds_; }

  template <class Node>
  bool CUSPARSECrsGraph<Node>::isInitialized() const
  { return isInitialized_; }

  template <class Scalar, class Node>
  CUSPARSECrsMatrix<Scalar,Node>::CUSPARSECrsMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<ParameterList> &params)
  : CrsMatrixBase<Scalar,int,Node>(graph,params) 
  , isInitialized_(false)
  {
    // Make sure that users only specialize for Kokkos Node types that are host Nodes (vs. device Nodes, such as GPU Nodes)
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Scalar, class Node>
  void CUSPARSECrsMatrix<Scalar,Node>::setValues(const ArrayRCP<const Scalar> &vals)
  { 
    std::string tfecfFuncName("setValues(vals)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true,
        std::runtime_error, " matrix is already initialized."
    )
    vals_ = vals;
    isInitialized_ = true;
  }

  template <class Scalar, class Node>
  ArrayRCP<const Scalar> CUSPARSECrsMatrix<Scalar,Node>::getValues() const
  { return vals_; }

  template <class Scalar, class Node>
  bool CUSPARSECrsMatrix<Scalar,Node>::isInitialized() const
  { return isInitialized_; }

  /// \class CUSPARSEOps
  /// \brief Default implementation of sparse matrix-vector multiply
  ///   and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  template <class Scalar, class Node>
  class CUSPARSEOps {
  public:
    //@{
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef int     ordinal_type;
    //! The Kokkos Node type.
    typedef Node    node_type;
    //! The type of this object, the sparse operator object
    typedef CUSPARSEOps<Scalar,Node> sparse_ops_type;

    /** \brief Typedef for local graph class */
    template <class O, class N>
    struct graph {
      typedef CUSPARSECrsGraph<N> graph_type;
    };

    /** \brief Typedef for local matrix class */
    template <class S, class O, class N>
    struct matrix {
      typedef CUSPARSECrsMatrix<S,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar type to the appropriate scalar.
    ///
    /// This always specifies a specialization of \c
    /// CUSPARSEOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef CUSPARSEOps<S2,Node> other_type;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    CUSPARSEOps(const RCP<Node> &node);

    //! Destructor
    ~CUSPARSEOps();

    //@}
    //! @name Accessor routines.
    //@{

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}

    //! @name Initialization of graph and matrix
    //@{

    //! \brief Allocate and initialize the storage for the matrix values.
    static ArrayRCP<size_t> allocRowPtrs(const ArrayView<const size_t> &rowPtrs);

    //! \brief Allocate and initialize the storage for a sparse graph.
    template <class T> 
    static ArrayRCP<T> allocStorage(const ArrayView<const size_t> &ptrs);

    //! Finalize a graph is null for CUSPARSE.
    static void finalizeGraph(CUSPARSECrsGraph<Node> &graph, const RCP<ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void finalizeMatrix(const CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params);
    
    //! Finalize a graph and a matrix.
    static void finalizeGraphAndMatrix(CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params);

    //! Initialize sparse operations with a graph and matrix
    void setGraphAndMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix);

    //@}
    //! @name Computational methods
    //@{

    /// \brief Y := alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, overwriting Y with the result.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result: \f$Y := \alpha A X\f$.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := Y + alpha * Op(A) * X.
    ///
    /// Apply the local sparse matrix A (or its transpose or conjugate
    /// transpose) to a multivector X, accumulating the result into Y.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to apply the matrix, its transpose,
    ///   or its conjugate transpose (if applicable).
    ///
    /// \param alpha [in] Scalar constant \f$\alpha\f$ by which to
    ///   multiply the result: \f$Y := Y + \alpha A X\f$.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [in/out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              RangeScalar beta,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Solve Y = Op(A) X for X, where we assume A is triangular.
    ///
    /// Solve the (upper or lower) triangular system Y = Op(A) X.
    /// Op(A) means A, the transpose of A, or the conjugate transpose
    /// of A, depending on the \c trans argument.
    ///
    /// \tparam DomainScalar The type of entries in the input
    ///   multivector X.  This may differ from the type of entries in
    ///   A or in Y.
    ///
    /// \tparam RangeScalar The type of entries in the output
    ///   multivector Y.  This may differ from the type of entries in
    ///   A or in X.
    ///
    /// \param trans [in] Whether to solve with the matrix, its
    ///   transpose, or its conjugate transpose (if applicable).
    ///
    /// \param uplo [in] UPPER_TRI if the matrix is upper triangular,
    ///   else LOWER_TRI if the matrix is lower triangular.
    ///
    /// \param diag [in] UNIT_DIAG if the matrix has an implicit unit diagonal,
    ///   else NON_UNIT_DIAG (diagonal entries are explicitly stored in the matrix).
    ///
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           Teuchos::EUplo uplo,
           Teuchos::EDiag diag,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    CUSPARSEOps(const CUSPARSEOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    RCP<Node> node_;

    size_t numRows_;
    bool isInitialized_;
           
    RCP<cusparseHandle_t> session_handle_; // same as Kokkos::CUSPARSEdetails::session_handle
    cusparseHybMat_t hybMat_;

    cusparseSolveAnalysisInfo_t ctransSolveAnalysis_, transSolveAnalysis_, notransSolveAnalysis_;
  };


  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::finalizeGraph(CUSPARSECrsGraph<Node> &graph, const RCP<ParameterList> &params)
  { 
    // nothing to do here
    std::string FuncName("Kokkos::CUSPARSEOps::finalizeGraph(graph,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false, 
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
  }

  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::finalizeMatrix(const CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params)
  { 
    // nothing much to do here
    std::string FuncName("Kokkos::CUSPARSEOps::finalizeMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false, 
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }

  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::finalizeGraphAndMatrix(CUSPARSECrsGraph<Node> &graph, CUSPARSECrsMatrix<Scalar,Node> &matrix, const RCP<ParameterList> &params)
  { 
    // nothing much to do here
    std::string FuncName("Kokkos::CUSPARSEOps::finalizeGraphAndMatrix(graph,matrix,params)");
    TEUCHOS_TEST_FOR_EXCEPTION(
        graph.isInitialized() == false, 
        std::runtime_error, FuncName << ": graph has not yet been initialized."
    )
    TEUCHOS_TEST_FOR_EXCEPTION(
        matrix.isInitialized() == false, 
        std::runtime_error, FuncName << ": matrix has not yet been initialized."
    )
  }


  template<class Scalar, class Node>
  CUSPARSEOps<Scalar,Node>::CUSPARSEOps(const RCP<Node> &node)
  : node_(node)
  , numRows_(0)
  , isInitialized_(false)
  {
    // Make sure that users only specialize CUSPARSEOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template<class Scalar, class Node>
  CUSPARSEOps<Scalar,Node>::~CUSPARSEOps() 
  {
    if (hybMat_) {
      cusparseDestroyHybMat(hybMat_);
    }
  }

  template <class Scalar, class Node>
  RCP<Node> CUSPARSEOps<Scalar,Node>::getNode() const {
    return node_;
  }

  template <class Scalar, class Node>
  void CUSPARSEOps<Scalar,Node>::setGraphAndMatrix(const RCP<const CUSPARSECrsGraph<Node> > &graph, const RCP<const CUSPARSECrsMatrix<Scalar,Node> > &matrix)
  {
    std::string tfecfFuncName("setGraphAndMatrix(graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        isInitialized_ == true, 
        std::runtime_error, " operators already initialized.");
    numRows_ = graph->getNumRows();
    // get cusparse objects from the matrix
    isInitialized_ = true;
  }

  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag,
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector< RangeScalar,Node> &X) const
  {
    std::string tfecfFuncName("solve(trans,uplo,diag,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPT(true)
    return;
  }


  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                      MultiVector<RangeScalar ,Node> &Y) const 
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,Y)");
    TEUCHOS_TEST_FOR_EXCEPT(true)
    return;
  }


  template <class Scalar, class Node>
  template <class DomainScalar, class RangeScalar>
  void CUSPARSEOps<Scalar,Node>::multiply(
                                Teuchos::ETransp trans,
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X,
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const
  {
    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPT(true)
    return;
  }


  // ======= pointer allocation ===========
  template <class Scalar, class Node>
  ArrayRCP<size_t> 
  CUSPARSEOps<Scalar,Node>::allocRowPtrs(const ArrayView<const size_t> &numEntriesPerRow) 
  {
    ArrayRCP<size_t> ptrs = arcp<size_t>( numEntriesPerRow.size() + 1 );
    ptrs[0] = 0;
    std::partial_sum( numEntriesPerRow.getRawPtr(), numEntriesPerRow.getRawPtr()+numEntriesPerRow.size(), ptrs.begin()+1 );
    return ptrs;
  }

  // ======= other allocation ===========
  template <class Scalar, class Node>
  template <class T>
  ArrayRCP<T> 
  CUSPARSEOps<Scalar,Node>::allocStorage(const ArrayView<const size_t> &rowPtrs)
  { 
    const size_t totalNumEntries = *(rowPtrs.end()-1);
    // alloc data
    ArrayRCP<T> vals;
    if (totalNumEntries > 0) vals = arcp<T>(totalNumEntries);
    std::fill( vals.begin(), vals.end(), Teuchos::ScalarTraits<T>::zero() );
    return vals;
  }

} // namespace Kokkos

#endif /* KOKKOS_CUSPARSEOPS_HPP */
