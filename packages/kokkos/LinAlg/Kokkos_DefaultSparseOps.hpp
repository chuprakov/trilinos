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

#ifndef KOKKOS_DEFAULTSPARSEOPS_HPP
#define KOKKOS_DEFAULTSPARSEOPS_HPP

#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_CompileTimeAssert.hpp>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrix.hpp" 
#include "Kokkos_CrsGraph.hpp" 
#include "Kokkos_MultiVector.hpp"
#include "Kokkos_NodeHelpers.hpp"
#include "Kokkos_DefaultArithmetic.hpp"
#include "Kokkos_DefaultSparseScaleKernelOps.hpp"
#include "Kokkos_DefaultSparseSolveKernelOps.hpp"
#include "Kokkos_DefaultSparseMultiplyKernelOps.hpp"

namespace Kokkos {

  /** \brief Default implementation of sparse matrix-vector multiplication and solve routines, for host-based nodes.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class DefaultHostSparseOps {
  public:
    //@{ 
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  ScalarType;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal OrdinalType;
    //! The Kokkos Node type.
    typedef Node    NodeType;

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The "rebind" struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// This always specifies a specialization of \c
    /// DefaultHostSparseOps, regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct rebind {
      typedef DefaultHostSparseOps<S2,Ordinal,Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    DefaultHostSparseOps(const RCP<Node> &node);

    //! Destructor
    ~DefaultHostSparseOps();

    //@}
    //! @name Accessor routines.
    //@{ 
    
    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Initialize structure of matrix, using CrsGraphHostCompute
    void initializeStructure(const CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &graph);

    //! Initialize values of matrix, using CrsMatrixHostCompute
    void initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;

    //! Left-scales a matrix by a vector
    template <class VectorScalar>
    void leftScale(const MultiVector<VectorScalar,Node> &X);

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultHostSparseOps(const DefaultHostSparseOps& source);

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> node_;

    // we do this one of two ways: 
    // 1D/packed: arrays of offsets, array of ordinals, array of values.
    ArrayRCP<const Ordinal> inds1D_;
    ArrayRCP<const size_t>  begs1D_, ends1D_;
    ArrayRCP<Scalar>  vals1D_;
    // 2D: array of pointers
    ArrayRCP<const ArrayRCP<Ordinal> > inds2D_;
    ArrayRCP<const ArrayRCP<Scalar> >  vals2D_;
    ArrayRCP<const size_t>          numEntries_;
    ArrayRCP<const Ordinal *> indPtrs_;
    ArrayRCP<Scalar  *> valPtrs_;

    size_t numRows_;
    bool indsInit_, valsInit_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultHostSparseOps<Scalar,Ordinal,Node>::DefaultHostSparseOps(const RCP<Node> &node)
  : node_(node) 
  {
    // Make sure that users only specialize DefaultHostSparseOps for
    // Kokkos Node types that are host Nodes (vs. device Nodes, such
    // as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultHostSparseOps<Scalar,Ordinal,Node>::~DefaultHostSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> DefaultHostSparseOps<Scalar,Ordinal,Node>::getNode() const {
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::clear() {
    begs1D_     = null;
    ends1D_     = null;
    inds1D_     = null;
    vals1D_     = null;
    inds2D_     = null;
    vals2D_     = null;
    numEntries_ = null;
    indPtrs_    = null;
    valPtrs_    = null;
    numRows_  = 0;
    indsInit_ = false;
    valsInit_ = false;
    isEmpty_  = false;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &graph) {
    using Teuchos::arcp;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else if (graph.is1DStructure()) {
      isEmpty_ = false;
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t> begs, ends;
      const_cast<CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(graph).get1DStructure( inds, begs, ends );
      inds1D_ = inds;
      begs1D_ = begs;
      ends1D_ = ends;
    }
    else {
      isEmpty_  = false;
      {
        ArrayRCP<ArrayRCP<Ordinal> > inds;
        ArrayRCP<size_t>  sizes;
        const_cast<CrsGraphHostCompute<Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(graph).get2DStructure(inds,sizes);
        inds2D_     = inds;
        numEntries_ = sizes;
      }
      indPtrs_    = arcp<const Ordinal *>(numRows_);
      for (size_t r=0; r < numRows_; ++r) {
        indPtrs_[r] = inds2D_[r].getRawPtr();
      }
    }
    indsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &matrix) {
    using Teuchos::arcp;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEUCHOS_TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty() || 
                       (inds2D_ != null && matrix.is1DStructure()) || (inds1D_ != null && matrix.is2DStructure()),
                       std::runtime_error, Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      if (matrix.is1DStructure()) {
        ArrayRCP<Scalar> vals;
        const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(matrix).get1DValues( vals );
        vals1D_ = vals;
      }
      else {
        {
          ArrayRCP<ArrayRCP<Scalar> > vals;
          const_cast<CrsMatrixHostCompute<Scalar,Ordinal,Node,DefaultHostSparseOps<void,Ordinal,Node> > &>(matrix).get2DValues(vals);
          vals2D_ = vals;
        }
        valPtrs_ = arcp<Scalar *>(numRows_);
        for (size_t r=0; r < numRows_; ++r) {
          valPtrs_[r] = vals2D_[r].getRawPtr();
        }
      }
    }
    valsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const 
  {
    typedef DefaultSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    typedef DefaultSparseSolveOp2<Scalar,Ordinal,DomainScalar,RangeScalar>  Op2D;
    typedef DefaultSparseTransposeSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp1D;
    typedef DefaultSparseTransposeSolveOp2<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp2D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): this solve was not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumCols() != Y.getNumCols(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumRows() < numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left-hand-side multivector does not have enough rows. Likely cause is that the column map was not provided to the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION(diag != Teuchos::UNIT_DIAG, std::runtime_error,
          Teuchos::typeName(*this) << "::solve(): solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else if (begs1D_ != null) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer< Scalar *>(valPtrs_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer< Scalar *>(valPtrs_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const 
  {
    // the 1 parameter to the template means that beta is ignored and the output multivector enjoys overwrite semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op1D;
    typedef DefaultSparseMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op2D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp1D;
    typedef DefaultSparseTransposeMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp2D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else if (begs1D_ != null) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        rbh.end();
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const 
  {
    // the 0 parameter means that beta is considered, and the output multivector enjoys accumulate semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op1D;
    typedef DefaultSparseMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op2D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp1D;
    typedef DefaultSparseTransposeMultiplyOp2<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp2D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y 
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else if (begs1D_ != null) {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
        wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op2D wdp;
        rbh.begin();
        wdp.numRows = numRows_;
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        rbh.end();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op2D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp2D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
        wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
        wdp.vals_beg   = rbh.template addConstBuffer<Scalar *>(valPtrs_);
        wdp.x          = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y          = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp2D>(0,numRHS,wdp);
      }
    }
    return;
  }




  template <class Scalar, class Ordinal, class Node>
  template <class VectorScalar>
  void DefaultHostSparseOps<Scalar,Ordinal,Node>::leftScale(const MultiVector<VectorScalar,Node> &X) 
  {
    typedef DefaultSparseScaleOp1<Scalar,Ordinal,VectorScalar>  Op1D;
    typedef DefaultSparseScaleOp2<Scalar,Ordinal,VectorScalar>  Op2D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::leftScale(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != 1);
    ReadyBufferHelper<Node> rbh(node_);
    if (begs1D_ != null) {
      Op1D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.begs    = rbh.template addConstBuffer<size_t>(begs1D_);
      wdp.ends    = rbh.template addConstBuffer<size_t>(ends1D_);
      wdp.inds    = rbh.template addConstBuffer<Ordinal>(inds1D_);
      wdp.vals    = rbh.template addNonConstBuffer<Scalar>(vals1D_);
      wdp.x       = rbh.template addConstBuffer<VectorScalar>(X.getValues());
      rbh.end();
      node_->template parallel_for<Op1D>(0,numRows_,wdp);	
    }
    else {
      Op2D wdp;
      rbh.begin();
      wdp.numRows = numRows_;
      wdp.numEntries = rbh.template addConstBuffer<size_t>(numEntries_);
      wdp.inds_beg   = rbh.template addConstBuffer<const Ordinal *>(indPtrs_);
      wdp.vals_beg   = rbh.template addNonConstBuffer<Scalar *>(valPtrs_);
      wdp.x          = rbh.template addConstBuffer<VectorScalar>(X.getValues());
      rbh.end();
      node_->template parallel_for<Op2D>(0,numRows_,wdp);
    }
    return;
  }


  /** \brief Default implementation of sparse matrix-vector multiplication and solve routines, for device-based nodes.
      \ingroup kokkos_crs_ops
    */
  template <class Scalar, class Ordinal, class Node>
  class DefaultDeviceSparseOps {
  public:
    //@{ 
    //! @name Typedefs and structs

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar  ScalarType;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal OrdinalType;
    //! The Kokkos Node type.
    typedef Node    NodeType;

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The "rebind" struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct rebind {
      typedef DefaultDeviceSparseOps<S2,Ordinal,Node> other;
    };

    //@}
    //! @name Constructors/Destructor
    //@{

    //! Constructor accepting and retaining a node object.
    DefaultDeviceSparseOps(const RCP<Node> &node);

    //! Destructor
    ~DefaultDeviceSparseOps();

    //@}
    //! @name Accessor routines.
    //@{ 
    
    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> getNode() const;

    //@}
    //! @name Initialization of structure
    //@{

    //! Initialize structure of matrix, using CrsGraphDeviceCompute
    void initializeStructure(const CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &graph);

    //! Initialize values of matrix, using CrsMatrixDeviceCompute
    void initializeValues(const CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &matrix);

    //! Clear all matrix structure and values.
    void clear();

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void multiply(Teuchos::ETransp trans, 
                  RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const;

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
    template <class DomainScalar, class RangeScalar>
    void solve(Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
               const MultiVector<DomainScalar,Node> &Y, MultiVector<RangeScalar,Node> &X) const;


    //! Left-scales a matrix by a vector
    template <class VectorScalar>
    void leftScale(const MultiVector<VectorScalar,Node> &X);

    //@}

  protected:
    //! Copy constructor (protected and unimplemented)
    DefaultDeviceSparseOps(const DefaultDeviceSparseOps& source);

    //! The Kokkos Node with which this object was instantiated.
    RCP<Node> node_;

    // we do this one of two ways: 
    // 1D: array of offsets, pointer for ordinals, pointer for values.
    ArrayRCP<const size_t>  pbuf_offsets1D_;
    ArrayRCP<const Ordinal> pbuf_inds1D_;
    ArrayRCP<Scalar>  pbuf_vals1D_;

    size_t numRows_;
    bool indsInit_, valsInit_, isEmpty_;
  };

  template<class Scalar, class Ordinal, class Node>
  DefaultDeviceSparseOps<Scalar,Ordinal,Node>::DefaultDeviceSparseOps(const RCP<Node> &node)
  : node_(node) 
  {
    clear();
  }

  template<class Scalar, class Ordinal, class Node>
  DefaultDeviceSparseOps<Scalar,Ordinal,Node>::~DefaultDeviceSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node>
  RCP<Node> DefaultDeviceSparseOps<Scalar,Ordinal,Node>::getNode() const { 
    return node_; 
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::clear() {
    pbuf_offsets1D_  = null;
    pbuf_inds1D_     = null;
    pbuf_vals1D_     = null;
    numRows_  = 0;
    indsInit_ = false;
    valsInit_ = false;
    isEmpty_  = false;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::initializeStructure(const CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &graph) {
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == true || valsInit_ == true, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeStructure(): structure already initialized.");
    numRows_ = graph.getNumRows();
    if (graph.isEmpty() || numRows_ == 0) {
      isEmpty_ = true;
    }
    else {
      ArrayRCP<Ordinal> inds;
      ArrayRCP<size_t > offs;
      const_cast<CrsGraphDeviceCompute<Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &>(graph).getDeviceBuffers(inds, offs);
      pbuf_inds1D_    = inds;
      pbuf_offsets1D_ = offs;
    }
    indsInit_ = true;
  }

  template <class Scalar, class Ordinal, class Node>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::initializeValues(const CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &matrix) {
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false, std::runtime_error, 
        Teuchos::typeName(*this) << "::initializeValues(): must initialize values after graph.");
    TEUCHOS_TEST_FOR_EXCEPTION(numRows_ != matrix.getNumRows() || isEmpty_ != matrix.isEmpty(), std::runtime_error,
        Teuchos::typeName(*this) << "::initializeValues(): matrix not compatible with previously supplied graph.");
    if (!isEmpty_) {
      ArrayRCP<Scalar> vals;
      const_cast<CrsMatrixDeviceCompute<Scalar,Ordinal,Node,DefaultDeviceSparseOps<void,Ordinal,Node> > &>(matrix).getDeviceBuffer(vals);
      pbuf_vals1D_ = vals;
    }
    valsInit_ = true;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::solve(
                      Teuchos::ETransp trans, Teuchos::EUplo uplo, Teuchos::EDiag diag, 
                      const MultiVector<DomainScalar,Node> &Y,
                            MultiVector<RangeScalar,Node> &X) const {
    typedef DefaultSparseSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  Op1D;
    typedef DefaultSparseTransposeSolveOp1<Scalar,Ordinal,DomainScalar,RangeScalar>  TOp1D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): this solve was not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumCols() != Y.getNumCols(), std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left hand side and right hand side multivectors have differing numbers of vectors.");
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumRows() < numRows_, std::runtime_error,
        Teuchos::typeName(*this) << "::solve(): Left-hand-side multivector does not have enough rows. Likely cause is that the column map was not provided to the Tpetra::CrsMatrix in the case of an implicit unit diagonal.");

    ReadyBufferHelper<Node> rbh(node_);
    if (numRows_ == 0) {
      // null op
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION(diag != Teuchos::UNIT_DIAG, std::runtime_error,
          Teuchos::typeName(*this) << "::solve(): solve of empty matrix only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign(X,Y);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addNonConstBuffer<DomainScalar>(X.getValuesNonConst());
        wdp.y       = rbh.template addConstBuffer<RangeScalar>(Y.getValues());
        rbh.end();
        wdp.numRows = numRows_;
        wdp.unitDiag = (diag == Teuchos::UNIT_DIAG ? true : false);
        wdp.upper    = (uplo == Teuchos::UPPER_TRI ? true : false);
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha,
                                const MultiVector<DomainScalar,Node> &X, 
                                MultiVector<RangeScalar,Node> &Y) const {
    // the 1 parameter to the template means that beta is ignored and the output multivector enjoys overwrite semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1>  Op1D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 1> TOp1D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= 0 * X
      //   <= 0
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Init(Y,Teuchos::ScalarTraits<RangeScalar>::zero());
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = Teuchos::ScalarTraits<RangeScalar>::zero(); // not used
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer<Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }


  template <class Scalar, class Ordinal, class Node>
  template <class DomainScalar, class RangeScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::multiply(
                                Teuchos::ETransp trans, 
                                RangeScalar alpha, const MultiVector<DomainScalar,Node> &X, 
                                RangeScalar beta, MultiVector<RangeScalar,Node> &Y) const {
    // the 0 parameter means that beta is considered, and the output multivector enjoys accumulate semantics
    typedef DefaultSparseMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0>  Op1D;
    typedef DefaultSparseTransposeMultiplyOp1<Scalar,Ordinal,DomainScalar,RangeScalar, 0> TOp1D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::multiply(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != Y.getNumCols());
    ReadyBufferHelper<Node> rbh(node_);
    if (isEmpty_ == true) {
      // Y <= alpha * 0 * X + beta * Y 
      //   <= beta * Y
      // NOTE: this neglects NaNs in X, which don't satisfy 0*NaN == 0
      //       however, X and Y may be of different size, and we need entries to determine how to mix those potential NaNs in X into Y
      //       therefore, the best we can do is scale Y to zero. Setting Y to zero would destroy NaNs in Y, which violates the semantics of the call.
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Scale(Y,beta);
    }
    else {
      if (trans == Teuchos::NO_TRANS) {
        Op1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
      }
      else {
        TOp1D wdp;
        rbh.begin();
        wdp.alpha   = alpha;
        wdp.beta    = beta;
        wdp.numRows = numRows_;
        wdp.numCols = Y.getNumRows();
        wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
        wdp.ends    = wdp.begs+1;
        wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
        wdp.vals    = rbh.template addConstBuffer< Scalar>(pbuf_vals1D_);
        wdp.x       = rbh.template addConstBuffer<DomainScalar>(X.getValues());
        wdp.y       = rbh.template addNonConstBuffer<RangeScalar>(Y.getValuesNonConst());
        wdp.xstride = X.getStride();
        wdp.ystride = Y.getStride();
        rbh.end();
        const size_t numRHS = X.getNumCols();
        node_->template parallel_for<TOp1D>(0,numRHS,wdp);
      }
    }
    return;
  }





  template <class Scalar, class Ordinal, class Node>
  template <class VectorScalar>
  void DefaultDeviceSparseOps<Scalar,Ordinal,Node>::leftScale(const MultiVector<VectorScalar,Node> &X)
  {
    typedef DefaultSparseScaleOp1<Scalar,Ordinal,VectorScalar>  Op1D;
    TEUCHOS_TEST_FOR_EXCEPTION(indsInit_ == false || valsInit_ == false, std::runtime_error,
        Teuchos::typeName(*this) << "::scale(): operation not fully initialized.");
    TEUCHOS_TEST_FOR_EXCEPT(X.getNumCols() != 1);
    ReadyBufferHelper<Node> rbh(node_);


    Op1D wdp;
    rbh.begin();
    wdp.numRows = numRows_;
    wdp.begs    = rbh.template addConstBuffer<size_t>(pbuf_offsets1D_);
    wdp.ends    = wdp.begs+1;
    wdp.inds    = rbh.template addConstBuffer<Ordinal>(pbuf_inds1D_);
    wdp.vals    = rbh.template addNonConstBuffer<Scalar>(pbuf_vals1D_);
    wdp.x       = rbh.template addConstBuffer<VectorScalar>(X.getValues());
    rbh.end();
    const size_t numRHS = X.getNumCols();
    node_->template parallel_for<Op1D>(0,numRows_*numRHS,wdp);
    
    return;
  }





  /** \example CrsMatrix_DefaultMultiplyTests.hpp 
    * This is an example that unit tests and demonstrates the implementation requirements for the DefaultSparseOps class.
    */

} // namespace Kokkos

#endif /* KOKKOS_DEFAULTSPARSEOPS_HPP */

