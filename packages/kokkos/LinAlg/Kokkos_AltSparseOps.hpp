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

#ifndef __Kokkos_AltSparseOps_hpp
#define __Kokkos_AltSparseOps_hpp

#include <Teuchos_CompileTimeAssert.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <iterator>
#include <stdexcept>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_CrsMatrixBase.hpp"
#include "Kokkos_CrsGraphBase.hpp"

#include "Kokkos_MultiVector.hpp"
#include "Kokkos_DefaultArithmetic.hpp"

#include "Kokkos_Raw_SparseMatVec_decl.hpp"
#include "Kokkos_Raw_SparseMatVec_def.hpp"
#include "Kokkos_Raw_SparseTriangularSolve_decl.hpp"
#include "Kokkos_Raw_SparseTriangularSolve_def.hpp"

#ifdef HAVE_KOKKOS_OPENMP
#  include "Kokkos_OpenMPNode.hpp"
#endif // HAVE_KOKKOS_OPENMP
#include "Kokkos_SerialNode.hpp"


namespace Kokkos {

  namespace details {

    /// \class AltSparseOpsDefaultAllocator
    /// \brief Default allocator for AltSparseOps.
    template <class Ordinal, class Node>
    class AltSparseOpsDefaultAllocator {
    public:
      static Teuchos::ArrayRCP<size_t>
      allocRowPtrs (const Teuchos::RCP<Node>& node,
                    const Teuchos::ArrayView<const size_t>& numEntriesPerRow)
      {
        (void) node;

        const Ordinal numRows = numEntriesPerRow.size ();

        // Let the ArrayRCP constructor do the allocation.
        Teuchos::ArrayRCP<size_t> ptr (numRows + 1);
        ptr[0] = 0; // We'll fill in the rest of the entries below.

        if (numRows > 0) {
          // Fill in ptr sequentially for now.  We might parallelize
          // this later, though it's still only O(numRows) work.
          // Parallel prefix is only likely to pay off for a large
          // number of threads.
          std::partial_sum (numEntriesPerRow.getRawPtr (),
                            numEntriesPerRow.getRawPtr () + numRows,
                            ptr.getRawPtr () + 1);
        }
        return ptr;
      }

      static Teuchos::ArrayRCP<size_t>
      copyRowPtrs (const Teuchos::RCP<Node>& node,
                   const Teuchos::ArrayView<const size_t>& rowPtrs)
      {
        (void) node;

        // Let the ArrayRCP constructor do the allocation.
        Teuchos::ArrayRCP<size_t> ptr (rowPtrs.size ());

        // Copy rowPtrs sequentially for now.  We might parallelize
        // this later, though it's still only O(numRows) work.
        std::copy (rowPtrs.getRawPtr (),
                   rowPtrs.getRawPtr () + rowPtrs.size (),
                   ptr.getRawPtr ());
        return ptr;
      }

      template<class T>
      static Teuchos::ArrayRCP<T>
      allocStorage (const Teuchos::RCP<Node>& node,
                    const Teuchos::ArrayView<const size_t>& rowPtrs)
      {
        (void) node;

        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsDefaultAllocator::allocStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const size_t totalNumEntries = rowPtrs[numRows];
        const T zero = Teuchos::ScalarTraits<T>::zero ();

        // Let the ArrayRCP constructor do the allocation.
        ArrayRCP<T> val (totalNumEntries);

        // Initialize the values sequentially.
        std::fill (val.getRawPtr (), val.getRawPtr () + totalNumEntries, zero);
        return val;
      }

      template<class T>
      static Teuchos::ArrayRCP<T>
      copyStorage (const Teuchos::RCP<Node>& node,
                   const Teuchos::ArrayView<const size_t>& rowPtrs,
                   const Teuchos::ArrayView<const T>& inputVals)
      {
        (void) node;

        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsDefaultAllocator::copyStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const size_t totalNumEntries = rowPtrs[numRows];

        TEUCHOS_TEST_FOR_EXCEPTION(
          inputVals.size() != totalNumEntries, std::invalid_argument,
          "AltSparseOpsDefaultAllocator::copyStorage: The inputVals array "
          "must have as many entries as the number of entries in the local "
          "sparse matrix, namely " << totalNumEntries << ".");

        // Let the ArrayRCP constructor do the allocation.
        ArrayRCP<T> val (totalNumEntries);

        // Get the raw pointers so that the compiler doesn't have to
        // optimize across the ArrayView inlining.
        const T* const rawInputVals = inputVals.getRawPtr ();
        T* const rawOutputVals = val.getRawPtr ();

        // Copy the values sequentially.
        std::copy (rawInputVals, rawInputVals + totalNumEntries, rawOutputVals);
        return val;
      }
    };


    // Kokkos Kernels for generic first-touch allocation.
    namespace {
      template<class Ordinal, class T>
      class ZeroInitKernel {
      private:
        T* const x_;

      public:
        ZeroInitKernel (T* const x) : x_ (x) {}

        void execute (const Ordinal i) {
          x_[i] = Teuchos::ScalarTraits<T>::zero ();
        }
      };

      template<class Ordinal, class T>
      class CopyKernel {
      private:
        T* const out_;
        const T* const in_;

      public:
        CopyKernel (T* const out, const T* const in) :
          out_ (out), in_ (in)
        {}

        void execute (const Ordinal i) {
          out_[i] = in_[i];
        }
      };

      template<class Ordinal, class T>
      class CsrInitKernel {
      private:
        const Ordinal* const ptr_;
        T* const val_;

      public:
        CsrInitKernel (const Ordinal* const ptr,
                       T* const val) :
          ptr_ (ptr),
          val_ (val)
        {}

        void execute (const Ordinal i) {
          for (Ordinal k = ptr_[i]; k < ptr_[i+1]; ++k) {
            val_[k] = Teuchos::ScalarTraits<T>::zero ();
          }
        }
      };


      template<class Ordinal, class T>
      class CsrCopyKernel {
      private:
        const Ordinal* const ptr_;
        T* const outVal_;
        const T* const inVal_;

      public:
        CsrCopyKernel (const Ordinal* const ptr,
                       T* const outVal,
                       const T* const inVal) :
          ptr_ (ptr),
          outVal_ (outVal),
          inVal_ (inVal)
        {}

        void execute (const Ordinal i) {
          for (Ordinal k = ptr_[i]; k < ptr_[i+1]; ++k) {
            outVal_[k] = inVal_[k];
          }
        }
      };
    } // namespace (anonymous)


    /// \class AltSparseOpsFirstTouchAllocator
    /// \brief Allocator for AltSparseOps, based on first-touch allocation.
    ///
    /// \tparam Ordinal The type for indexing into (local) data.
    /// \tparam Node The Kokkos Node type.
    template <class Ordinal, class Node>
    class AltSparseOpsFirstTouchAllocator {
    public:
      static Teuchos::ArrayRCP<size_t>
      allocRowPtrs (const Teuchos::RCP<Node>& node,
                    const Teuchos::ArrayView<const size_t>& numEntriesPerRow)
      {
        using Teuchos::ArrayRCP;
        using Teuchos::arcp;
        const Ordinal numRows = numEntriesPerRow.size ();

        // Allocate raw, since arcp() might initialize in debug mode.
        size_t* const rawPtr = new size_t [numRows + 1];

        // Initialize the row pointers in parallel, using the Kokkos
        // Node.  If the Kokkos Node's parallelization scheme respects
        // first-touch initialization, this should set the proper NUMA
        // affinity, at least at page boundaries.
        typedef ZeroInitKernel<Ordinal, Ordinal> kernel_type;
        node->parallel_for (0, numRows+1, kernel_type (rawPtr));

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        ArrayRCP<Ordinal> ptr = arcp<Ordinal> (rawPtr, 0, numRows+1, true);

        if (numRows > 0) {
          // Fill in ptr sequentially for now.  We might parallelize
          // this later, though it's still only O(numRows) work.
          // Parallel prefix is only likely to pay off for a large
          // number of threads.
          std::partial_sum (numEntriesPerRow.getRawPtr (),
                            numEntriesPerRow.getRawPtr () + numRows,
                            ptr.getRawPtr () + 1);
        }
        return ptr;
      }

      static Teuchos::ArrayRCP<Ordinal>
      copyRowPtrs (const Teuchos::RCP<Node>& node,
                   const Teuchos::ArrayView<const Ordinal>& rowPtrs)
      {
        using Teuchos::arcp;
        const Ordinal numRows = rowPtrs.size () - 1;

        // Don't force the compiler to inline ArrayView::operator[] in
        // the loop below.
        const Ordinal* const rawRowPtrs = rowPtrs.getRawPtr ();

        // Allocate raw, since arcp() might initialize in debug mode.
        Ordinal* const rawPtr = new Ordinal [numRows + 1];

        // Copy the row pointers in parallel.  If first touch works,
        // this should set the proper affinity at page boundaries.
        typedef CopyKernel<Ordinal, Ordinal> kernel_type;
        node->parallel_for (0, numRows+1, kernel_type (rawPtr, rawRowPtrs));

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<Ordinal> (rawPtr, 0, numRows+1, true);
      }

      template<class T>
      static Teuchos::ArrayRCP<T>
      allocStorage (const Teuchos::RCP<Node>& node,
                    const Teuchos::ArrayView<const Ordinal>& rowPtrs)
      {
        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::allocStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const Ordinal totalNumEntries = rowPtrs[numRows];
        const T zero = Teuchos::ScalarTraits<T>::zero ();

        // Allocate raw, since arcp() might initialize in debug mode.
        T* const rawVal = new T [totalNumEntries];

        // Get the row pointer so that the compiler doesn't have to
        // optimize the ArrayView.
        const Ordinal* const rawRowPtrs = rowPtrs.getRawPtr ();

        // Initialize the values in parallel.  If first touch works,
        // this should set the proper affinity at page boundaries.  We
        // iterate over the values using the row pointers, so that if
        // the Kokkos Node parallelizes in a reproducible way, it will
        // initialize the values using the same affinity with which
        // the row pointers were initialized.
        CsrInitKernel<Ordinal, T> kernel;
        node->parallel_for (0, numRows, kernel);

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<T> (rawVal, 0, totalNumEntries, true);
      }


      template<class T>
      static Teuchos::ArrayRCP<T>
      copyStorage (const Teuchos::RCP<Node>& node,
                   const Teuchos::ArrayView<const Ordinal>& rowPtrs,
                   const Teuchos::ArrayView<const T>& inputVals)
      {
        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::copyStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const Ordinal totalNumEntries = rowPtrs[numRows];

        TEUCHOS_TEST_FOR_EXCEPTION(
          inputVals.size() != totalNumEntries, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::copyStorage: The inputVals array "
          "must have as many entries as the number of entries in the local "
          "sparse matrix, namely " << totalNumEntries << ".");

        // Allocate raw, since arcp() might initialize in debug mode.
        T* const rawOutputVals = new T [totalNumEntries];

        // Get the raw pointers so that the compiler doesn't have to
        // optimize across the ArrayView inlining.
        const Ordinal* const rawRowPtrs = rowPtrs.getRawPtr ();
        const T* const rawInputVals = inputVals.getRawPtr ();

        // Copy the values in parallel.  If first touch works, this
        // should set the proper affinity at page boundaries.  We
        // iterate over the values using the row pointers, so that if
        // the Kokkos Node parallelizes in a reproducible way, it will
        // initialize the values using the same affinity with which
        // the row pointers were initialized.
        typedef CsrCopyKernel<Ordinal, T> kernel_type;
        kernel_type kernel (rawRowPtrs, rawOutputVals, rawInputVals);
        node->parallel_for (0, numRows, kernel);

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<T> (rawOutputVals, 0, totalNumEntries, true);
      }
    };


#ifdef HAVE_KOKKOS_OPENMP
    // Partial speicalization of AltSparseOpsFirstTouchAllocator for
    // Node=Kokkos::OpenMPNode.  This just uses OpenMP pragmas
    // directly, avoiding any possible overhead of going through the
    // Kokkos Node instance.
    template <class Ordinal>
    class AltSparseOpsFirstTouchAllocator<Ordinal, OpenMPNode> {
    public:
      static Teuchos::ArrayRCP<Ordinal>
      allocRowPtrs (const Teuchos::RCP<OpenMPNode>& node,
                    const Teuchos::ArrayView<const Ordinal>& numEntriesPerRow)
      {
        (void) node;

        using Teuchos::ArrayRCP;
        using Teuchos::arcp;
        const Ordinal numRows = numEntriesPerRow.size ();

        // Allocate raw, since arcp() might initialize in debug mode.
        Ordinal* const rawPtr = new Ordinal [numRows + 1];

        // Initialize the row pointers in parallel.  If the operating
        // system respects first-touch initialization, this should set
        // the proper NUMA affinity, at least at page boundaries.  We
        // don't rely on the Kokkos Node for parallelization; instead,
        // we use OpenMP, since that's what our kernels use.  If
        // you're building without OpenMP support, the compiler will
        // simply ignore this pragma.
        #pragma omp parallel for
        for (Ordinal i = 0; i < numRows+1; ++i) {
          rawPtr[i] = 0;
        }
        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        ArrayRCP<size_t> ptr = arcp<size_t> (rawPtr, 0, numRows+1, true);

        if (numRows > 0) {
          // Fill in ptr sequentially for now.  We might parallelize
          // this later, though it's still only O(numRows) work.
          // Parallel prefix is only likely to pay off for a large
          // number of threads.
          std::partial_sum (numEntriesPerRow.getRawPtr (),
                            numEntriesPerRow.getRawPtr () + numRows,
                            ptr.getRawPtr () + 1);
        }
        return ptr;
      }

      static Teuchos::ArrayRCP<Ordinal>
      copyRowPtrs (const Teuchos::RCP<OpenMPNode>& node,
                   const Teuchos::ArrayView<const size_t>& rowPtrs)
      {
        (void) node;

        using Teuchos::arcp;
        const Ordinal numRows = rowPtrs.size () - 1;

        // Don't force the compiler to inline ArrayView::operator[] in
        // the loop below.
        const size_t* const rawRowPtrs = rowPtrs.getRawPtr ();

        // Allocate raw, since arcp() might initialize in debug mode.
        size_t* const rawPtr = new size_t [numRows + 1];

        // Copy the row pointers in parallel.  If first touch works,
        // this should set the proper affinity at page boundaries.  We
        // don't rely on the Kokkos Node for this; instead, we use
        // OpenMP, since that's what our kernels use.  If you're
        // building without OpenMP support, the compiler will simply
        // ignore this pragma.
        #pragma omp parallel for
        for (Ordinal i = 0; i < numRows+1; ++i) {
          rawPtr[i] = rawRowPtrs[i];
        }

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<size_t> (rawPtr, 0, numRows+1, true);
      }

      template<class T>
      static Teuchos::ArrayRCP<T>
      allocStorage (const Teuchos::RCP<OpenMPNode>& node,
                    const Teuchos::ArrayView<const size_t>& rowPtrs)
      {
        (void) node;

        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::allocStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const size_t totalNumEntries = rowPtrs[numRows];
        const T zero = Teuchos::ScalarTraits<T>::zero ();

        // Allocate raw, since arcp() might initialize in debug mode.
        T* const rawVal = new T [totalNumEntries];

        // Get the row pointer so that the compiler doesn't have to
        // optimize the ArrayView.
        const size_t* const rawRowPtrs = rowPtrs.getRawPtr ();

        // Initialize the values in parallel.  If first touch works,
        // this should set the proper affinity at page boundaries.  We
        // don't rely on the Kokkos Node for this; instead, we use
        // OpenMP, since that's what our kernels use.  We also iterate
        // over the values using the row pointers, so that if OpenMP
        // parallelizes in a reproducible way, it will initialize the
        // values using the same affinity with which the row pointers
        // were initialized.
        //
        // If you're building without OpenMP support, the compiler will
        // simply ignore this pragma.
#pragma omp parallel for
        for (Ordinal i = 0; i < numRows; ++i) {
          for (size_t k = rawRowPtrs[i]; k < rawRowPtrs[i+1]; ++k) {
            rawVal[k] = zero;
          }
        }

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<T> (rawVal, 0, totalNumEntries, true);
      }


      template<class T>
      static Teuchos::ArrayRCP<T>
      copyStorage (const Teuchos::RCP<Kokkos::OpenMPNode>& node,
                   const Teuchos::ArrayView<const size_t>& rowPtrs,
                   const Teuchos::ArrayView<const T>& inputVals)
      {
        (void) node;

        using Teuchos::ArrayRCP;
        using Teuchos::arcp;

        TEUCHOS_TEST_FOR_EXCEPTION(rowPtrs.size() == 0, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::copyStorage: The input rowPtrs array "
          "must have length at least one, but rowPtrs.size() = 0.  (Remember that "
          "rowPtrs must have exactly one more entry than the number of rows in "
          "the matrix.)");

        const Ordinal numRows = rowPtrs.size() - 1;
        const size_t totalNumEntries = rowPtrs[numRows];

        TEUCHOS_TEST_FOR_EXCEPTION(
          inputVals.size() != totalNumEntries, std::invalid_argument,
          "AltSparseOpsFirstTouchAllocator::copyStorage: The inputVals array "
          "must have as many entries as the number of entries in the local "
          "sparse matrix, namely " << totalNumEntries << ".");

        // Allocate raw, since arcp() might initialize in debug mode.
        T* const rawOutputVals = new T [totalNumEntries];

        // Get the raw pointers so that the compiler doesn't have to
        // optimize across the ArrayView inlining.
        const size_t* const rawRowPtrs = rowPtrs.getRawPtr ();
        const T* const rawInputVals = inputVals.getRawPtr ();

        // Copy the values in parallel.  If first touch works, this
        // should set the proper affinity at page boundaries.  We
        // don't rely on the Kokkos Node for this; instead, we use
        // OpenMP, since that's what our kernels use.  We also iterate
        // over the values using the row pointers, so that if OpenMP
        // parallelizes in a reproducible way, it will initialize the
        // values using the same affinity with which the row pointers
        // were initialized.
        //
        // If you're building without OpenMP support, the compiler will
        // simply ignore this pragma.
#pragma omp parallel for
        for (Ordinal i = 0; i < numRows; ++i) {
          for (size_t k = rawRowPtrs[i]; k < rawRowPtrs[i+1]; ++k) {
            rawOutputVals[k] = rawInputVals[k];
          }
        }

        // Encapsulate the raw pointer in an (owning) ArrayRCP.
        return arcp<T> (rawOutputVals, 0, totalNumEntries, true);
      }
    };
#endif // HAVE_KOKKOS_OPENMP

  } // namespace details

  /// \class AltCrsGraph
  /// \brief Local sparse graph in compressed sparse row format;
  ///   suitable for host-based Kokkos Nodes.
  template <class Ordinal,class Node>
  class AltCrsGraph {
  public:
    typedef Ordinal ordinal_type;
    typedef Node node_type;

    AltCrsGraph (Ordinal numRows, Ordinal numCols,
                 const Teuchos::RCP<Node>& node,
                 const Teuchos::RCP<Teuchos::ParameterList>& params);

    Teuchos::RCP<Node> getNode() const {
      return node_;
    }

    Ordinal getNumRows () const {
      return numRows_;
    }

    Ordinal getNumCols () const {
      return numCols_;
    }

    Teuchos::ArrayRCP<const size_t> getPointers() const {
      return ptr_;
    }

    Teuchos::ArrayRCP<const Ordinal> getIndices() const {
      return ind_;
    }

    bool isInitialized() const {
      return isInitialized_;
    }

    /// \brief Whether the graph is empty.
    ///
    /// "Empty" means either that the graph has no rows (the number of
    /// rows is zero), or that the graph has no stored entries.
    bool isEmpty() const {
      return isEmpty_;
    }

    /// \brief Whether the graph has empty rows.
    ///
    /// An empty graph (see isEmpty()) trivially has empty rows.
    /// Otherwise, the graph has empty rows if one or more rows
    /// contains no stored entries.
    bool hasEmptyRows() const {
      return hasEmptyRows_;
    }

    void setStructure (const Teuchos::ArrayRCP<const size_t>& ptr,
                       const Teuchos::ArrayRCP<const Ordinal>& ind);
    void setMatDesc (Teuchos::EUplo uplo, Teuchos::EDiag diag);
    void getMatDesc (Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const;

  private:
    Teuchos::RCP<Node> node_;
    Ordinal numRows_, numCols_;
    //Teuchos::RCP<ParameterList> params_;
    bool isInitialized_;
    bool isEmpty_;
    bool hasEmptyRows_;
    Teuchos::EUplo tri_uplo_;
    Teuchos::EDiag unit_diag_;

    Teuchos::ArrayRCP<const size_t> ptr_;
    Teuchos::ArrayRCP<const Ordinal> ind_;
  };


  /// \class AltCrsMatrix
  /// \brief Local sparse matrix in compressed sparse row format;
  ///   suitable for host-based Kokkos Nodes.
  ///
  /// \note Tied to a particular AltCrsGraph instance that defines the
  ///   structure of the sparse matrix.
  template <class Scalar,
            class Ordinal,
            class Node>
  class AltCrsMatrix {
  public:
    typedef Scalar scalar_type;
    typedef Ordinal ordinal_type;
    typedef Node node_type;
    typedef AltCrsGraph<Ordinal,Node> graph_type;

    AltCrsMatrix (const Teuchos::RCP<const AltCrsGraph<Ordinal,Node> > &graph,
                  const Teuchos::RCP<Teuchos::ParameterList> &params);

    void setValues (const Teuchos::ArrayRCP<const Scalar>& val);

    Teuchos::ArrayRCP<const Scalar> getValues() const {
      return val_;
    }

    bool isInitialized() const {
      return isInitialized_;
    }

  private:
    Teuchos::RCP<const graph_type> graph_;
    Teuchos::ArrayRCP<const Scalar> val_;
    bool isInitialized_;
  };

  template <class Ordinal, class Node>
  AltCrsGraph<Ordinal,Node>::
  AltCrsGraph (Ordinal numRows, Ordinal numCols,
               const Teuchos::RCP<Node> &node,
               const Teuchos::RCP<Teuchos::ParameterList> &params) :
    node_ (node),
    numRows_ (numRows),
    numCols_ (numCols),
    isInitialized_ (false),
    isEmpty_ (numRows == 0), // provisional; a matrix with numRows > 0
                             // may still have zero entries.
    hasEmptyRows_ (true), // provisional
    tri_uplo_ (Teuchos::UNDEF_TRI),
    unit_diag_ (Teuchos::NON_UNIT_DIAG)
  {
    // Make sure that users only specialize for Kokkos Node types that
    // are host Nodes (vs. device Nodes, such as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void) cta;

    // We don't use params currently.
    (void) params;
  }

  template <class Ordinal, class Node>
  void
  AltCrsGraph<Ordinal,Node>::
  setStructure (const Teuchos::ArrayRCP<const size_t> &ptr,
                const Teuchos::ArrayRCP<const Ordinal> &ind)
  {
    std::string tfecfFuncName("setStructure(ptr,ind)");
    const Ordinal numRows = this->getNumRows();

    // mfh 19 June 2012: The tests expect std::runtime_error rather
    // than the arguably more appropriate std::invalid_argument, so
    // I'll throw std::runtime_error here.  Ditto for the other checks
    // below.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr.is_null (),
      std::runtime_error,
      ": The input array 'ptr' must be nonnull, even for a matrix with zero "
      "rows.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr.size() != (size_t) numRows+1,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptr.size() = " << ptr.size() << " != numRows+1 = "
      << (numRows+1) << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ptr[0] != 0,
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ptr[0] = " << ptr[0] << " != 0.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ind.size() != (size_t) ptr[numRows],
      std::runtime_error,
      ": Graph input data are not coherent:\n"
      "-- ind.size() = " << ind.size() << " != ptr[numRows="
      << numRows << "] = " << ptr[numRows] << ".");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_,
      std::runtime_error,
      ": Graph has already been initialized."
    )

    const size_t numEntries = ptr[numRows];
    if (numRows == 0 || numEntries == 0) {
      isEmpty_ = true;
      hasEmptyRows_ = true; // trivially
    }
    else {
      isEmpty_ = false;
      // Check whether the graph has any empty rows.
      bool emptyRows = false;
      for (Ordinal i = 0; i < numRows; ++i) {
        if (ptr[i] == ptr[i+1]) {
          emptyRows = true;
          break;
        }
      }
      hasEmptyRows_ = emptyRows;
    }
    ptr_ = ptr;
    ind_ = ind;
    isInitialized_ = true;
  }

  template <class Ordinal, class Node>
  void
  AltCrsGraph<Ordinal,Node>::
  setMatDesc (Teuchos::EUplo uplo, Teuchos::EDiag diag)
  {
    tri_uplo_ = uplo;
    unit_diag_ = diag;
  }

  template <class Ordinal, class Node>
  void
  AltCrsGraph<Ordinal,Node>::
  getMatDesc (Teuchos::EUplo &uplo, Teuchos::EDiag &diag) const
  {
    uplo = tri_uplo_;
    diag = unit_diag_;
  }

  template <class Scalar, class Ordinal, class Node>
  AltCrsMatrix<Scalar,Ordinal,Node>::
  AltCrsMatrix (const Teuchos::RCP<const AltCrsGraph<Ordinal,Node> >& graph,
                const Teuchos::RCP<Teuchos::ParameterList>& params) :
    graph_ (graph),
    isInitialized_ (false)
  {
    // Make sure that users only specialize for Kokkos Node types that
    // are host Nodes (vs. device Nodes, such as GPU Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void) cta;

    // We don't use params currently.
    (void) params;
  }

  template <class Scalar, class Ordinal, class Node>
  void
  AltCrsMatrix<Scalar,Ordinal,Node>::
  setValues (const Teuchos::ArrayRCP<const Scalar> &val)
  {
    std::string tfecfFuncName("setValues(val)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_, std::runtime_error, ": The matrix is already initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      graph_.is_null() && ! val.is_null(),
      std::runtime_error,
      ": The matrix has a null graph, but you're trying to give it a nonnull "
      "array of values."
    );
    val_ = val;
    if (val_.is_null ()) {
      isInitialized_ = false;
    }
    else {
      isInitialized_ = true;
    }
  }

  /// \class AltSparseOps
  /// \brief Alternate implementation of local sparse matrix-vector
  ///   multiply and solve routines, for host-based Kokkos Node types.
  /// \ingroup kokkos_crs_ops
  ///
  /// \tparam Scalar The type of entries of the sparse matrix.
  /// \tparam Ordinal The type of (local) indices of the sparse matrix.
  /// \tparam Node The Kokkos Node type.
  /// \tparam Allocator Class that defines static methods for
  ///   allocating the various arrays used by the compressed sparse
  ///   row format.  This is where first-touch allocation can be
  ///   implemented, for example.
  ///
  /// This class is called AltSparseOps ("alternate sparse
  /// operations") because the default local sparse operations class
  /// on host-based Kokkos Nodes is DefaultHostSparseOps.
  ///
  /// While this class is templated on the Kokkos Node type, it does
  /// not use any functionality of the Kokkos Node.
  ///
  /// This class uses the given Allocator by default.
  template <class Scalar,
            class Ordinal,
            class Node,
            class Allocator = details::AltSparseOpsDefaultAllocator<Ordinal, Node> >
  class AltSparseOps : public Teuchos::Describable {
  public:
    //! \name Typedefs and structs
    //@{

    //! The type of the individual entries of the sparse matrix.
    typedef Scalar scalar_type;
    //! The type of the (local) indices describing the structure of the sparse matrix.
    typedef Ordinal ordinal_type;
    //! The Kokkos Node type.
    typedef Node node_type;
    //! The Allocator type, whose class methods allocate the three arrays.
    typedef Allocator allocator_type;
    //! The type of this object, the sparse operator object.
    typedef AltSparseOps<Scalar, Ordinal, Node, Allocator> sparse_ops_type;

    //! Typedef for local graph class
    template <class O, class N>
    struct graph {
      typedef AltCrsGraph<O,N> graph_type;
    };

    //! Typedef for local matrix class
    template <class S, class O, class N>
    struct matrix {
      typedef AltCrsMatrix<S,O,N> matrix_type;
    };

    /// \brief Sparse operations type for a different scalar type.
    ///
    /// The bind_scalar struct defines the type responsible for sparse
    /// operations for a scalar type S2, which may be different from
    /// \c Scalar.
    ///
    /// Use by Tpetra CrsMatrix to bind a potentially "void" scalar
    /// type to the appropriate scalar.
    ///
    /// This always specifies a specialization of AltSparseOps,
    /// regardless of the scalar type S2.
    ///
    /// \tparam S2 A scalar type possibly different from \c Scalar.
    template <class S2>
    struct bind_scalar {
      typedef AltSparseOps<S2, Ordinal, Node> other_type;
    };

    //@}
    //! \name Constructors and destructor
    //@{

    /// \brief Constructor, with default parameters.
    ///
    /// We syntactically forbid setting parameters after construction,
    /// since setting parameters after calling setGraphAndMatrix()
    /// would require reorganizing the already optimized sparse matrix
    /// storage.  If you want to set nondefault values of parameters,
    /// you must use the constructor that takes a ParameterList.
    ///
    /// \param node [in/out] Kokkos Node instance.
    AltSparseOps (const Teuchos::RCP<Node>& node);

    /// \brief Constructor, with custom parameters.
    ///
    /// Both this constructor and finalizeGraphAndMatrix() accept a
    /// ParameterList.  However, those sets of parameters are
    /// different.  The constructor's parameters concern the
    /// algorithm, and the parameters for finalizeGraphAndMatrix()
    /// concern the data structure.  It's possible to use different
    /// algorithms with the same data structure.
    ///
    /// \param node [in/out] Kokkos Node instance.
    ///
    /// \param params [in/out] Parameters for the solve.  We fill in
    ///   the given ParameterList with its default values, but we
    ///   don't keep it around.  (This saves a bit of memory.)
    AltSparseOps (const Teuchos::RCP<Node>& node,
                  Teuchos::ParameterList& plist);

    //! Destructor
    ~AltSparseOps();

    /// \brief Get a default ParameterList.
    ///
    /// The returned ParameterList has all accepted parameters, their
    /// default values, documentation, and validators (if applicable).
    ///
    /// This is a class (static) method so that you can get the
    /// default ParameterList (with built-in documentation) before
    /// constructing a AltSparseOps instance.
    static Teuchos::RCP<const Teuchos::ParameterList>
    getDefaultParameters ()
    {
      using Teuchos::ParameterList;
      using Teuchos::parameterList;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;

      RCP<ParameterList> plist = parameterList ("AltSparseOps");
      setDefaultParameters (*plist);
      return rcp_const_cast<const ParameterList> (plist);
    }

    //@}
    //! \name Implementation of Teuchos::Describable
    //@{

    //! One-line description of this instance.
    std::string description () const {
      using Teuchos::TypeNameTraits;
      std::ostringstream os;
      os <<  "Kokkos::AltSparseOps<"
         << "Scalar=" << TypeNameTraits<Scalar>::name()
         << ", Ordinal=" << TypeNameTraits<Ordinal>::name()
         << ", Node=" << TypeNameTraits<Node>::name()
         << ">";
      return os.str();
    }

    //! Write a possibly more verbose description of this instance to out.
    void
    describe (Teuchos::FancyOStream& out,
              const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const
    {
      using Teuchos::EVerbosityLevel;
      using Teuchos::includesVerbLevel;
      using Teuchos::OSTab;
      using Teuchos::rcpFromRef;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;
      using std::endl;

      // Interpret the default verbosity level as VERB_MEDIUM.
      const EVerbosityLevel vl =
        (verbLevel == VERB_DEFAULT) ? VERB_MEDIUM : verbLevel;

      if (vl == VERB_NONE) {
        return;
      }
      else if (includesVerbLevel (vl, VERB_LOW)) { // vl >= VERB_LOW
        out << this->description();

        if (includesVerbLevel (vl, VERB_MEDIUM)) { // vl >= VERB_MEDIUM
          out << ":" << endl;
          OSTab tab1 (rcpFromRef (out));

          out << "matVecVariant_ = " << matVecVariant_ << endl
              << "unroll_ = " << unroll_ << endl
              << "isInitialized_ = " << isInitialized_ << endl;
          if (isInitialized_) {
            std::string triUplo ("INVALID");
            if (tri_uplo_ == Teuchos::UNDEF_TRI) {
              triUplo = "UNDEF_TRI";
            }
            else if (tri_uplo_ == Teuchos::LOWER_TRI) {
              triUplo = "LOWER_TRI";
            }
            else if (tri_uplo_ == Teuchos::UPPER_TRI) {
              triUplo = "UPPER_TRI";
            }
            std::string unitDiag ("INVALID");
            if (unit_diag_ == Teuchos::NON_UNIT_DIAG) {
              unitDiag = "NON_UNIT_DIAG";
            }
            else if (unit_diag_ == Teuchos::UNIT_DIAG) {
              unitDiag = "UNIT_DIAG";
            }

            out << "numRows_ = " << numRows_ << endl
                << "isEmpty_ = " << isEmpty_ << endl
                << "hasEmptyRows_ = " << hasEmptyRows_ << endl
                << "tri_uplo_ = " << triUplo << endl
                << "unit_diag_ = " << unitDiag << endl;
            if (ptr_.size() > 0) {
              out << "numEntries = " << ptr_[ptr_.size()-1] << endl;
            }
            else {
              out << "numEntries = 0" << endl;
            }

            if (includesVerbLevel (vl, VERB_EXTREME)) { // vl >= VERB_EXTREME
              // Only print out all the sparse matrix's data in
              // extreme verbosity mode.
              out << "ptr_ = [";
              std::copy (ptr_.begin(), ptr_.end(),
                         std::ostream_iterator<Ordinal> (out, " "));
              out << "]" << endl << "ind_ = [";
              std::copy (ind_.begin(), ind_.end(),
                         std::ostream_iterator<Ordinal> (out, " "));
              out << "]" << endl << "val_ = [";
              std::copy (val_.begin(), val_.end(),
                         std::ostream_iterator<Scalar> (out, " "));
              out << "]" << endl;
            } // vl >= VERB_EXTREME
          } // if is initialized
        } // vl >= VERB_MEDIUM
      } // vl >= VERB_LOW
    }

    //@}
    //! \name Accessor routines
    //@{

    //! The Kokkos Node with which this object was instantiated.
    Teuchos::RCP<Node> getNode () const;

    //@}
    //! @name Initialization of graph and matrix
    //@{

    /// \brief Allocate and initialize the storage for the row pointers.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param numEntriesPerRow [in] numEntriesPerRow[i] is the number
    ///   of entries in row i, for all rows of the local sparse matrix.
    static Teuchos::ArrayRCP<size_t>
    allocRowPtrs (const Teuchos::RCP<Node>& node,
                  const ArrayView<const size_t>& numEntriesPerRow);

    /// \brief Copy the storage for the row pointers.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param rowPtrs [in] The array of row offsets to copy.
    ///
    /// You might like to call this method if the Allocator promises a
    /// special allocation method (like first-touch allocation), but
    /// the input array was not allocated by the Allocator.
    static Teuchos::ArrayRCP<Ordinal>
    copyRowPtrs (const Teuchos::RCP<Node>& node,
                 const Teuchos::ArrayView<const Ordinal>& rowPtrs);

    /// \brief Allocate and initialize the storage for a sparse graph or matrix.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param rowPtrs [in] The array of row offsets; the 'ptr' array
    ///   in the compressed sparse row storage format.  rowPtrs.size()
    ///   is one plus the number of rows in the local sparse matrix.
    template <class T>
    static Teuchos::ArrayRCP<T>
    allocStorage (const Teuchos::RCP<Node>& node,
                  const Teuchos::ArrayView<const size_t>& rowPtrs);

    /// \brief Copy the storage for a sparse graph or matrix.
    ///
    /// \param node [in/out] Kokkos Node instance.
    /// \param rowPtrs [in] The array of row offsets; the 'ptr' array
    ///   in the compressed sparse row storage format.  rowPtrs.size()
    ///   is one plus the number of rows in the local sparse matrix.
    /// \param inputVals [in] The array of input values (or column
    ///   indices) to copy.
    ///
    /// You might like to call this method if the Allocator promises a
    /// special allocation method (like first-touch allocation), but
    /// inputVals was not allocated by the Allocator.
    template <class T>
    static Teuchos::ArrayRCP<T>
    copyStorage (const Teuchos::RCP<Node>& node,
                 const Teuchos::ArrayView<const Ordinal>& rowPtrs,
                 const Teuchos::ArrayView<const T>& inputVals);

    //! Finalize the graph.
    static void
    finalizeGraph (Teuchos::EUplo uplo,
                   Teuchos::EDiag diag,
                   AltCrsGraph<Ordinal, Node>& graph,
                   const Teuchos::RCP<Teuchos::ParameterList> &params);

    //! Finalize the matrix of an already-finalized graph.
    static void
    finalizeMatrix (const AltCrsGraph<Ordinal, Node>& graph,
                    AltCrsMatrix<Scalar, Ordinal, Node>& matrix,
                    const Teuchos::RCP<Teuchos::ParameterList>& params);

    /// \brief Finalize a graph and a matrix.
    ///
    /// Both the constructor and this method accept a ParameterList.
    /// However, those sets of parameters are different.  The
    /// constructor's parameters concern the algorithm, and the
    /// parameters for this method concern the data structure.  It's
    /// possible to use different algorithms with the same data
    /// structure.
    static void
    finalizeGraphAndMatrix (Teuchos::EUplo uplo,
                            Teuchos::EDiag diag,
                            AltCrsGraph<Ordinal, Node>& graph,
                            AltCrsMatrix<Scalar, Ordinal, Node>& matrix,
                            const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Initialize sparse operations with a graph and matrix.
    void
    setGraphAndMatrix (const Teuchos::RCP<const AltCrsGraph<Ordinal,Node> > &graph,
                       const Teuchos::RCP<const AltCrsMatrix<Scalar,Ordinal,Node> > &matrix);

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
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param Y [out] Result multivector. Contents will be overwritten.
    template <class DomainScalar, class RangeScalar>
    void
    multiply (Teuchos::ETransp trans,
              RangeScalar alpha,
              const MultiVector<DomainScalar,Node> &X,
              MultiVector<RangeScalar,Node> &Y) const;

    /// \brief Y := beta * Y + alpha * Op(A) * X.
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
    ///   multiply the result of the sparse matrix-(multi)vector
    ///   multiply.
    ///
    /// \param X [in] Input multivector.
    ///
    /// \param beta [in] Scalar constant \f$\beta\f$ by which to
    ///   multiply Y when summing with the result of the sparse
    ///   matrix-(multi)vector multiply.
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
    /// \param Y [in] Input multivector.
    ///
    /// \param X [out] Result multivector.
    template <class DomainScalar, class RangeScalar>
    void
    solve (Teuchos::ETransp trans,
           const MultiVector<DomainScalar,Node> &Y,
           MultiVector<RangeScalar,Node> &X) const;

    //@}

  private:
    /// \brief Which algorithm variant to use for sparse matrix-vector multiply.
    ///
    /// The textbook compressed sparse row (CSR) and compressed sparse
    /// column (CSC) sparse matrix-vector multiply algorithms have two
    /// nested 'for' loops.  The outer for loop is for the rows (for
    /// CSR; columns for CSC), and the inner for loop is for the
    /// entries within a row (for CSR; column for CSC).  We call this
    /// the 'for-for' variant.
    ///
    /// We also make available two variants that use a single 'for'
    /// loop over all the entries in the sparse matrix.  The first,
    /// which we call 'for-while', has an inner whlie loop for
    /// incrementing the current row (for CSR; column for CSC) index.
    /// The second, which we call 'for-if', replaces the while loop in
    /// 'for-while' with a single if statement.  The 'for-if' variant
    /// is only correct if the sparse matrix contains no empty rows
    /// (for CSR; columns for CSC).  If you specify the for-if
    /// variant, we check for empty rows, and use the for-while
    /// variant if there are any empty rows.
    enum EMatVecVariant {
      FOR_FOR,
      FOR_WHILE,
      FOR_IF
    };

    //! Fill the given ParameterList with defaults and validators.
    static void
    setDefaultParameters (Teuchos::ParameterList& plist);

    /// \brief Set the default mat-vec algorithm variant parameter.
    ///
    /// Use this to construct a ParameterList with default values and
    /// validators.
    static void
    setDefaultMatVecVariantParameter (Teuchos::ParameterList& plist);

    //! Set the default multivector unroll parameter.
    static void
    setDefaultUnrollParameter (Teuchos::ParameterList& plist);

    //! Set the default parameter for forcing first-touch allocation.
    static void
    setDefaultFirstTouchParameter (Teuchos::ParameterList& plist);

    //! Copy constructor (protected and unimplemented)
    AltSparseOps (const AltSparseOps& source);

    //! The Kokkos Node instance given to this object's constructor.
    Teuchos::RCP<Node> node_;

    //! \name Raw compressed sparse row storage of the local sparse matrix.
    //@{
    ArrayRCP<const size_t>  ptr_;
    ArrayRCP<const Ordinal> ind_;
    ArrayRCP<const Scalar>  val_;
    //@}

    Teuchos::EUplo  tri_uplo_;
    Teuchos::EDiag unit_diag_;

    //! \name Parameters set by the constructor's input ParameterList.
    //@{
    EMatVecVariant matVecVariant_;
    bool unroll_;
    bool firstTouchAllocation_;
    //@}

    Ordinal numRows_;
    bool isInitialized_;
    bool isEmpty_;
    bool hasEmptyRows_;
  };

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  finalizeGraph (Teuchos::EUplo uplo,
                 Teuchos::EDiag diag,
                 AltCrsGraph<Ordinal,Node>& graph,
                 const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isInitialized(), std::runtime_error,
      "Kokkos::AltSparseOps::finalizeGraph: Graph has not yet been initialized."
    );
    graph.setMatDesc (uplo, diag);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  finalizeMatrix (const AltCrsGraph<Ordinal,Node> &graph,
                  AltCrsMatrix<Scalar,Ordinal,Node> &matrix,
                  const Teuchos::RCP<Teuchos::ParameterList> &params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! matrix.isInitialized(), std::runtime_error,
      "Kokkos::AltSparseOps::finalizeMatrix(graph,matrix,params): "
      "Matrix has not yet been initialized."
    );
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  finalizeGraphAndMatrix (Teuchos::EUplo uplo,
                          Teuchos::EDiag diag,
                          AltCrsGraph<Ordinal,Node>& graph,
                          AltCrsMatrix<Scalar,Ordinal,Node>& matrix,
                          const Teuchos::RCP<Teuchos::ParameterList>& params)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! graph.isInitialized(), std::runtime_error,
      "Kokkos::AltSparseOps::finalizeGraphAndMatrix(graph,matrix,params): "
      "Graph has not yet been initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! matrix.isInitialized(), std::runtime_error,
      "Kokkos::AltSparseOps::finalizeGraphAndMatrix(graph,matrix,params): "
      "Matrix has not yet been initialized."
    );
    graph.setMatDesc (uplo, diag);
  }


  template <class Scalar, class Ordinal, class Node, class Allocator>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  AltSparseOps (const Teuchos::RCP<Node>& node) :
    node_ (node),
    tri_uplo_ (Teuchos::UNDEF_TRI),      // Provisionally
    unit_diag_ (Teuchos::NON_UNIT_DIAG), // Provisionally
    matVecVariant_ (AltSparseOps<Scalar, Ordinal, Node, Allocator>::FOR_FOR),
    unroll_ (true),
    firstTouchAllocation_ (false),
    numRows_ (0),                        // Provisionally
    isInitialized_ (false),
    isEmpty_ (true),                     // Provisionally
    hasEmptyRows_ (true)                 // Provisionally
  {
    // Make sure that users only specialize AltSparseOps for Kokkos
    // Node types that are host Nodes (vs. device Nodes, such as GPU
    // Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  AltSparseOps (const Teuchos::RCP<Node>& node,
                Teuchos::ParameterList& params) :
    node_ (node),
    tri_uplo_ (Teuchos::UNDEF_TRI),      // Provisionally
    unit_diag_ (Teuchos::NON_UNIT_DIAG), // Provisionally
    matVecVariant_ (AltSparseOps<Scalar, Ordinal, Node, Allocator>::FOR_FOR),
    unroll_ (true),
    firstTouchAllocation_ (false),
    numRows_ (0),                        // Provisionally
    isInitialized_ (false),
    isEmpty_ (true),                     // Provisionally
    hasEmptyRows_ (true)                 // Provisionally
  {
    using Teuchos::getIntegralValue;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // Make sure that users only specialize AltSparseOps for Kokkos
    // Node types that are host Nodes (vs. device Nodes, such as GPU
    // Nodes).
    Teuchos::CompileTimeAssert<Node::isHostNode == false> cta; (void)cta;

    params.validateParametersAndSetDefaults (*getDefaultParameters ());

    const std::string varParamName ("Sparse matrix-vector multiply variant");
    const std::string unrollParamName ("Unroll across multivectors");
    const std::string firstTouchParamName ("Force first-touch allocation");

    matVecVariant_ = getIntegralValue<EMatVecVariant> (params, varParamName);
    unroll_ = params.get<bool> (unrollParamName);
    firstTouchAllocation_ = params.get<bool> (firstTouchParamName);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::~AltSparseOps() {
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  setDefaultParameters (Teuchos::ParameterList& plist)
  {
    setDefaultMatVecVariantParameter (plist);
    setDefaultUnrollParameter (plist);
    setDefaultFirstTouchParameter (plist);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  setDefaultUnrollParameter (Teuchos::ParameterList& plist)
  {
    const bool unroll = true;
    plist.set ("Unroll across multivectors", unroll, "Whether to unroll reads "
               "and writes of multivectors across columns of the input and "
               "ouput multivectors");
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  setDefaultFirstTouchParameter (Teuchos::ParameterList& plist)
  {
    const bool firstTouchAllocation = false;
    plist.set ("Force first-touch allocation", firstTouchAllocation,
               "Whether to copy all the input data and initialize them in "
               "parallel (using OpenMP if available), to ensure first-touch "
               "initialization");
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  setDefaultMatVecVariantParameter (Teuchos::ParameterList& plist)
  {
    using Teuchos::Array;
    using Teuchos::ParameterList;
    using Teuchos::stringToIntegralParameterEntryValidator;

    Array<std::string> strs (3);
    strs[0] = "for-for";
    strs[1] = "for-while";
    strs[2] = "for-if";

    Array<std::string> docs (3);
    docs[0] = "Two nested for loops (textbook algorithm)";
    docs[1] = "Outer for loop, inner while loop";
    docs[2] = "Outer for loop, inner if statement";

    Array<EMatVecVariant> vals (3);
    vals[0] = FOR_FOR;
    vals[1] = FOR_WHILE;
    vals[2] = FOR_IF;

    const std::string paramName ("Sparse matrix-vector multiply variant");
    const std::string paramDoc ("Which algorithm variant to use for sparse "
                                "matrix-vector multiply");
    plist.set (paramName, strs[0], paramDoc,
               stringToIntegralParameterEntryValidator<EMatVecVariant> (strs(), docs(), vals(), strs[0]));
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  RCP<Node> AltSparseOps<Scalar, Ordinal, Node, Allocator>::getNode() const {
    return node_;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  setGraphAndMatrix (const Teuchos::RCP<const AltCrsGraph<Ordinal,Node> > &opgraph,
                     const Teuchos::RCP<const AltCrsMatrix<Scalar,Ordinal,Node> > &opmatrix)
  {
    using Teuchos::ArrayRCP;
    using Teuchos::arcp;
    using Teuchos::arcp_const_cast;
    using Teuchos::null;
    // using std::cerr;
    // using std::endl;

    std::string tfecfFuncName("setGraphAndMatrix(uplo,diag,graph,matrix)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isInitialized_, std::runtime_error, " operators already initialized."
    );

    ArrayRCP<const  size_t> ptr = opgraph->getPointers ();
    ArrayRCP<const Ordinal> ind = opgraph->getIndices ();
    ArrayRCP<const Scalar> val = opmatrix->getValues ();
    const Ordinal numRows = opgraph->getNumRows ();

    // Verify the input data before setting internal state.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr.size() != (size_t) numRows + 1,
      std::invalid_argument,
      ": ptr.size() = " << ptr.size() << " != numRows+1 = " << (numRows + 1) << ".");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ind.size() != val.size(),
      std::invalid_argument,
      ": ind.size() = " << ind.size() << " != val.size() = " << val.size()
      << ", for ptr = opgraph->getPointers() and ind = opgraph->getIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      (size_t) ptr[numRows] != (size_t) ind.size(),
      std::invalid_argument,
      ": ptr[numRows = " << numRows << "] = " << ptr[numRows]
      << " != ind.size() = " << ind.size() << ", for ptr = "
      "opgraph->getPointers() and ind = opgraph->getIndices().");

    numRows_ = numRows;
    hasEmptyRows_ = opgraph->hasEmptyRows ();

    if (opgraph->isEmpty () || numRows_ == 0) {
      isEmpty_ = true;
      // We have to go through a little trouble because ptr_ is an
      // array of const size_t, but we need to set its first entry
      // here.
      ArrayRCP<size_t> myPtr = arcp<size_t> (1);
      myPtr[0] = 0;
      ptr_ = arcp_const_cast<const size_t> (myPtr);
      myPtr = null;

      // The matrix is empty, so ind_ and val_ have zero entries.
      ind_ = null;
      val_ = null;
    }
    else {
      isEmpty_ = false;

      if (firstTouchAllocation_) {
        // Override the input data's allocation.  This assumes that
        // the input graph and matrix were not allocated using
        // first-touch allocation.  Make copies, initializing them
        // using the first-touch procedure.  This also overrides the
        // Allocator template parameter.
        typedef details::AltSparseOpsFirstTouchAllocator<Ordinal, Node> alloc_type;

        ptr_ = alloc_type::copyRowPtrs (node_, ptr ());
        ind_ = alloc_type::template copyStorage<Ordinal> (node_, ptr (), ind ());
        val_ = alloc_type::template copyStorage<Scalar> (node_, ptr (), val ());
      }
      else { // No first touch allocation.
        // We can just use these arrays directly.
        ptr_ = ptr;
        ind_ = ind;
        val_ = val;
      }
    }
    opgraph->getMatDesc (tri_uplo_, unit_diag_);
    isInitialized_ = true;
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  solve (Teuchos::ETransp trans,
         const MultiVector<DomainScalar, Node>& Y,
         MultiVector<RangeScalar, Node>& X) const
  {
    std::string tfecfFuncName("solve(trans,Y,X)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isInitialized_,
      std::runtime_error,
      ": The solve was not fully initialized."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumCols() != Y.getNumCols(),
      std::runtime_error,
      ": Input and output multivectors have different numbers of vectors (columns)."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      static_cast<size_t> (X.getNumRows()) < static_cast<size_t> (numRows_),
      std::runtime_error,
      ": Output multivector X does not have enough rows.  X.getNumRows() == "
      << X.getNumRows() << ", but the matrix has " << numRows_ << " rows."
    );

    if (numRows_ == 0) {
      return; // Nothing to do
    }
    else if (isEmpty_) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        unit_diag_ != Teuchos::UNIT_DIAG,
        std::runtime_error,
        ": Solve with empty matrix is only valid for an implicit unit diagonal.");
      // solve I * X = Y for X = Y
      DefaultArithmetic<MultiVector<RangeScalar,Node> >::Assign (X, Y);
    }
    else {
      typedef ordinal_type OT;
      typedef scalar_type MST; // matrix scalar type
      typedef DomainScalar DST;
      typedef RangeScalar RST;

      RST* const X_raw = X.getValuesNonConst ().getRawPtr ();
      const Ordinal X_stride = (Ordinal) X.getStride ();
      const DST* const Y_raw = Y.getValues ().getRawPtr ();
      const Ordinal Y_stride = (Ordinal) Y.getStride ();

      const  size_t* const ptr = ptr_.getRawPtr ();
      const Ordinal* const ind = ind_.getRawPtr ();
      const Scalar*  const val = val_.getRawPtr ();
      const Ordinal numRows = X.getNumRows ();
      const Ordinal numCols = Y.getNumRows ();
      const Ordinal numVecs = X.getNumCols ();

      if (trans == Teuchos::NO_TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCsrColMajorUnitDiag;
            lowerTriSolveCsrColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else { // non unit diagonal
            using Kokkos::Raw::lowerTriSolveCsrColMajor;
            lowerTriSolveCsrColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCsrColMajorUnitDiag;
            upperTriSolveCsrColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else { // non unit diagonal
            using Kokkos::Raw::upperTriSolveCsrColMajor;
            upperTriSolveCsrColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCscColMajorUnitDiag;
            // numRows resp. numCols come from the number of rows in Y
            // resp. X, so they still appear in the same order as
            // in the not transposed cases above.
            lowerTriSolveCscColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::lowerTriSolveCscColMajor;
            lowerTriSolveCscColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCscColMajorUnitDiag;
            upperTriSolveCscColMajorUnitDiag<OT, MST, DST, RST> (numRows, numCols,
                                                                 numVecs,
                                                                 X_raw, X_stride,
                                                                 ptr, ind, val,
                                                                 Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::upperTriSolveCscColMajor;
            upperTriSolveCscColMajor<OT, MST, DST, RST> (numRows, numCols,
                                                         numVecs,
                                                         X_raw, X_stride,
                                                         ptr, ind, val,
                                                         Y_raw, Y_stride);
          }
        }
      }
      else if (trans == Teuchos::CONJ_TRANS) {
        if (tri_uplo_ == Teuchos::LOWER_TRI) {
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::lowerTriSolveCscColMajorUnitDiagConj;
            lowerTriSolveCscColMajorUnitDiagConj<OT, MST, DST, RST> (numRows, numCols,
                                                                     numVecs,
                                                                     X_raw, X_stride,
                                                                     ptr, ind, val,
                                                                     Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::lowerTriSolveCscColMajorConj;
            lowerTriSolveCscColMajorConj<OT, MST, DST, RST> (numRows, numCols,
                                                             numVecs,
                                                             X_raw, X_stride,
                                                             ptr, ind, val,
                                                             Y_raw, Y_stride);
          }
        }
        else { // upper triangular
          if (unit_diag_ == Teuchos::UNIT_DIAG) {
            using Kokkos::Raw::upperTriSolveCscColMajorUnitDiagConj;
            upperTriSolveCscColMajorUnitDiagConj<OT, MST, DST, RST> (numRows, numCols,
                                                                     numVecs,
                                                                     X_raw, X_stride,
                                                                     ptr, ind, val,
                                                                     Y_raw, Y_stride);
          }
          else {
            using Kokkos::Raw::upperTriSolveCscColMajorConj;
            upperTriSolveCscColMajorConj<OT, MST, DST, RST> (numRows, numCols,
                                                             numVecs,
                                                             X_raw, X_stride,
                                                             ptr, ind, val,
                                                             Y_raw, Y_stride);
          }
        }
      }
    }
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar, Node>& X,
            MultiVector<RangeScalar, Node>& Y) const
  {
    typedef DomainScalar DST;
    typedef RangeScalar RST;
    const RST beta = Teuchos::ScalarTraits<RST>::zero ();
    this->template multiply<DST, RST> (trans, alpha, X, beta, Y);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class DomainScalar, class RangeScalar>
  void
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  multiply (Teuchos::ETransp trans,
            RangeScalar alpha,
            const MultiVector<DomainScalar, Node> &X,
            RangeScalar beta,
            MultiVector<RangeScalar, Node> &Y) const
  {
    using std::cerr;
    using std::endl;

    std::string tfecfFuncName("multiply(trans,alpha,X,beta,Y)");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isInitialized_,
      std::runtime_error,
      ": Sparse ops not initialized.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      X.getNumCols() != Y.getNumCols(),
      std::runtime_error,
      ": X and Y do not have the same number of columns.");

    typedef ordinal_type OT;
    typedef scalar_type MST; // matrix scalar type
    typedef DomainScalar DST;
    typedef RangeScalar RST;

    // These dimensions come from the input and output multivectors,
    // so they apply for the transpose case as well.
    const Ordinal numRows = Y.getNumRows ();
    const Ordinal numCols = X.getNumRows ();
    const Ordinal numVecs = X.getNumCols ();
    RST* const Y_raw = Y.getValuesNonConst ().getRawPtr ();
    const Ordinal Y_stride = Y.getStride ();
    const DST* const X_raw = X.getValues ().getRawPtr ();
    const Ordinal X_stride = X.getStride ();

    const  size_t* const ptr = ptr_.getRawPtr ();
    const Ordinal* const ind = ind_.getRawPtr ();
    const Scalar*  const val = val_.getRawPtr ();

    // Pointer to the sparse matrix-vector multiply routine to use.
    void (*matVec) (const OT, const OT, const OT,
                    const RST&, RST[], const OT,
                    const RST&, const size_t[], const OT[], const MST[],
                    const DST[], const OT);

    // The following very long switch statement selects one of the
    // hard-coded-numVecs routines for certain values of numVecs.
    // (Hard-coding the number of columns in the multivectors avoids
    // two branches and an integer addition.)  Otherwise, it picks a
    // general routine.  Here is also where we use the parameters
    // given to the constructor to pick the algorithm variant.
    //
    // Note that we're taking numRows and numCols from Y resp. X.
    // Assuming that the dimensions of X and Y are correct, then
    // whether or not we're applying the transpose, the (transposed,
    // if applicable) matrix has dimensions numRows by numCols.
    // That's why we don't switch the order of numRows, numCols in the
    // invocations below.
    switch (numVecs) {
    case 1:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor1Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif1Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile1Vec<OT, MST, DST, RST>;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor1Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif1Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile1Vec<OT, MST, DST, RST>;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj1Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj1Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj1Vec<OT, MST, DST, RST>;
        }
      }
      break;
    case 2:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor2Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif2Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile2Vec<OT, MST, DST, RST>;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor2Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif2Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile2Vec<OT, MST, DST, RST>;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj2Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj2Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj2Vec<OT, MST, DST, RST>;
        }
      }
      break;
    case 3:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor3Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif3Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile3Vec<OT, MST, DST, RST>;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor3Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif3Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile3Vec<OT, MST, DST, RST>;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj3Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj3Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj3Vec<OT, MST, DST, RST>;
        }
      }
      break;
    case 4:
      if (trans == Teuchos::NO_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForfor4Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCsrColMajorForif4Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCsrColMajorForwhile4Vec<OT, MST, DST, RST>;
        }
      }
      else if (trans == Teuchos::TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForfor4Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForif4Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhile4Vec<OT, MST, DST, RST>;
        }
      }
      else { // if (trans == Teuchos::CONJ_TRANS) {
        if (matVecVariant_ == FOR_FOR) {
          matVec = &Kokkos::Raw::matVecCscColMajorForforConj4Vec<OT, MST, DST, RST>;
        }
        else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
          matVec = &Kokkos::Raw::matVecCscColMajorForifConj4Vec<OT, MST, DST, RST>;
        }
        else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
          matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj4Vec<OT, MST, DST, RST>;
        }
      }
      break;
    default: // The "general case"
      if (unroll_) {
        if (trans == Teuchos::NO_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForfor4Unrolled<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            // FIXME (mfh 26 Jul 2012) Currently, the code generator
            // is broken for the 'for-if' and 'for-while' variants,
            // when the number of multivector columns is not fixed (it
            // wasn't broken before -- it had to do with the
            // introduction of temporaries).  We could either revert
            // to the for-for variant, or throw an exception.  We do
            // the latter, since for-if is not the default, and we
            // want to make sure that people don't get confusing
            // performance results.

            //matVec = &Kokkos::Raw::matVecCsrColMajorForif4Unrolled<OT, MST, DST, RST>;
            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE ||
                 // (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCsrColMajorForwhile4Unrolled<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-while' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
        else if (trans == Teuchos::TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForfor4Unrolled<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            //matVec = &Kokkos::Raw::matVecCscColMajorForif4Unrolled<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE || (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCscColMajorForwhile4Unrolled<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-while' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
        else { // if (trans == Teuchos::CONJ_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForforConj4Unrolled<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            //matVec = &Kokkos::Raw::matVecCscColMajorForifConj4Unrolled<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE ||
                 // (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj4Unrolled<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-while' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
      }
      else { // Don't unroll across multivector columns
        if (trans == Teuchos::NO_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCsrColMajorForfor<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            //matVec = &Kokkos::Raw::matVecCsrColMajorForif<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE ||
                 // (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCsrColMajorForwhile<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
        else if (trans == Teuchos::TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForfor<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            //matVec = &Kokkos::Raw::matVecCscColMajorForif<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE ||
                 // (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCscColMajorForwhile<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-while' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
        else { // if (trans == Teuchos::CONJ_TRANS) {
          if (matVecVariant_ == FOR_FOR) {
            matVec = &Kokkos::Raw::matVecCscColMajorForforConj<OT, MST, DST, RST>;
          }
          else if (matVecVariant_ == FOR_IF && ! hasEmptyRows_) {
            //matVec = &Kokkos::Raw::matVecCscColMajorForifConj<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-if' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
          else { // matVecVariant_ == FOR_WHILE ||
                 // (matVecVariant_ == FOR_IF && hasEmptyRows_)
            //matVec = &Kokkos::Raw::matVecCscColMajorForwhileConj<OT, MST, DST, RST>;

            TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "The 'for-while' "
              "variant of sparse matrix-vector multiply is not currently "
              "implemented for a non-fixed number of columns in the multi"
              "vectors.  Please use the 'for-for' variant for now.");
          }
        }
      }
    }

    // Now we know what mat-vec routine to call, so call it.
    matVec (numRows, numCols, numVecs, beta, Y_raw, Y_stride,
            alpha, ptr, ind, val, X_raw, X_stride);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  Teuchos::ArrayRCP<size_t>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  allocRowPtrs (const Teuchos::RCP<Node>& node,
                const Teuchos::ArrayView<const size_t>& numEntriesPerRow)
  {
    return allocator_type::allocRowPtrs (node, numEntriesPerRow);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  Teuchos::ArrayRCP<Ordinal>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  copyRowPtrs (const Teuchos::RCP<Node>& node,
               const Teuchos::ArrayView<const Ordinal>& rowPtrs)
  {
    return allocator_type::copyRowPtrs (node, rowPtrs);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class T>
  Teuchos::ArrayRCP<T>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  allocStorage (const Teuchos::RCP<Node>& node,
                const Teuchos::ArrayView<const size_t>& rowPtrs)
  {
    return allocator_type::template allocStorage<T> (node, rowPtrs);
  }

  template <class Scalar, class Ordinal, class Node, class Allocator>
  template <class T>
  Teuchos::ArrayRCP<T>
  AltSparseOps<Scalar, Ordinal, Node, Allocator>::
  copyStorage (const Teuchos::RCP<Node>& node,
               const Teuchos::ArrayView<const Ordinal>& rowPtrs,
               const Teuchos::ArrayView<const T>& inputVals)
  {
    return allocator_type::template allocStorage<T> (node, rowPtrs, inputVals);
  }

} // namespace Kokkos

#endif // #ifndef __Kokkos_AltSparseOps_hpp

