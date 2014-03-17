// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

// FINISH: can't use rowPtrs_ without checking that it exists
// FINISH: add code to fillComplete() and CrsMatrix::fillComplete() to delete the Tpetra data

#ifndef TPETRA_KOKKOSREFACTOR_CRSGRAPH_DEF_HPP
#define TPETRA_KOKKOSREFACTOR_CRSGRAPH_DEF_HPP

#include <Tpetra_Distributor.hpp>
#include <Teuchos_Assert.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_NullIteratorTraits.hpp>
#include <Teuchos_as.hpp>
#include <algorithm>
#include <string>
#include <utility>
#include <Teuchos_SerialDenseMatrix.hpp>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_CrsGraph_decl.hpp"
#endif
#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

namespace Tpetra {

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            size_t maxNumEntriesPerRow,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_ (true)
  , noRedundancies_ (true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumEntriesPerRow == OTST::invalid(),
      std::invalid_argument, "The allocation hint must be a valid size_t value, "
      "which in this case means it must not be Teuchos::OrdinalTraits<size_t>::"
      "invalid().");
    resumeFill(params);
    checkInternalState();
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &colMap,
            size_t maxNumEntriesPerRow,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocForAllRows_(maxNumEntriesPerRow)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_(true)
  , noRedundancies_(true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION(maxNumEntriesPerRow == OTST::invalid(),
      std::invalid_argument, "The allocation hint must be a valid size_t value, "
      "which in this case means it must not be Teuchos::OrdinalTraits<size_t>::"
      "invalid().");
    resumeFill(params);
    checkInternalState();
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_(true)
  , noRedundancies_(true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    const char tfecfFuncName[] = "CrsGraph(rowMap,NumEntriesPerRowToAlloc)";
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::invalid_argument,
        ": NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    for (size_t r=0; r < getNodeNumRows(); ++r) {
      const size_t curRowCount = NumEntriesPerRowToAlloc[r];
      TEUCHOS_TEST_FOR_EXCEPTION(curRowCount == OTST::invalid(),
        std::invalid_argument, "NumEntriesPerRowToAlloc[" << r << "] specifies "
        "an invalid number of entries (Teuchos::OrdinalTraits<size_t>::"
        "invalid()).");
    }
    resumeFill(params);
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &colMap,
            const ArrayRCP<const size_t> &NumEntriesPerRowToAlloc,
            ProfileType pftype,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(pftype)
  , numAllocPerRow_(NumEntriesPerRowToAlloc)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(false)
  , indicesAreLocal_(false)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_(true)
  , noRedundancies_(true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    typedef Teuchos::OrdinalTraits<size_t> OTST;
    const char tfecfFuncName[] = "CrsGraph(rowMap,colMap,NumEntriesPerRowToAlloc)";
    staticAssertions();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)NumEntriesPerRowToAlloc.size() != getNodeNumRows(), std::invalid_argument,
        ": NumEntriesPerRowToAlloc must have as many entries as specified by rowMap for this node.");
    for (size_t r=0; r < getNodeNumRows(); ++r) {
      const size_t curRowCount = NumEntriesPerRowToAlloc[r];
      TEUCHOS_TEST_FOR_EXCEPTION(curRowCount == OTST::invalid(),
        std::invalid_argument, "NumEntriesPerRowToAlloc[" << r << "] specifies "
        "an invalid number of entries (Teuchos::OrdinalTraits<size_t>::"
        "invalid()).");
    }
    resumeFill(params);
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &colMap,
            const t_RowPtrs & rowPointers,
            const t_LocalOrdinal_1D & columnIndices,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(StaticProfile)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(true)
  , indicesAreLocal_(true)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_(true)
  , noRedundancies_(true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    staticAssertions();
    globalNumEntries_ = globalNumDiags_ = globalMaxNumRowEntries_ = OrdinalTraits<global_size_t>::invalid();
    setAllIndices(rowPointers,columnIndices);
    checkInternalState();
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rowMap,
            const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &colMap,
            const ArrayRCP<size_t> & rowPointers,
            const ArrayRCP<LocalOrdinal> & columnIndices,
            const RCP<ParameterList>& params)
  : DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >(rowMap)
  , rowMap_(rowMap)
  , colMap_(colMap)
  , nodeNumEntries_(0)
  , nodeNumAllocated_(OrdinalTraits<size_t>::invalid())
  , pftype_(StaticProfile)
  , numAllocForAllRows_(0)
  , indicesAreAllocated_(true)
  , indicesAreLocal_(true)
  , indicesAreGlobal_(false)
  , fillComplete_(false)
  , indicesAreSorted_(true)
  , noRedundancies_(true)
  , haveLocalConstants_ (false)
  , haveGlobalConstants_ (false)
  , haveRowInfo_(true)
  {
    staticAssertions();
    globalNumEntries_ = globalNumDiags_ = globalMaxNumRowEntries_ = OrdinalTraits<global_size_t>::invalid();
    setAllIndices(rowPointers,columnIndices);
    checkInternalState();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
    typename KokkosClassic::DefaultKernels<
      void,
      LocalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  CrsGraph (const RCP<const map_type>& rowMap,
            const RCP<const map_type>& colMap,
            const LocalStaticCrsGraphType& k_local_graph_,
            const RCP<ParameterList>& params)
    : DistObject<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, node_type> (rowMap)
    , rowMap_ (rowMap)
    , colMap_ (colMap)
    , k_lclGraph_ (k_local_graph_)
    , globalNumEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalNumDiags_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , globalMaxNumRowEntries_ (Teuchos::OrdinalTraits<global_size_t>::invalid ())
    , nodeNumEntries_ (0) // FIXME (mfh 17 Mar 2014) should get from k_lclGraph_ right now
    , nodeNumAllocated_ (Teuchos::OrdinalTraits<size_t>::invalid ())
    , pftype_ (StaticProfile)
    , numAllocForAllRows_ (0)
    , indicesAreAllocated_ (true)
    , indicesAreLocal_ (true)
    , indicesAreGlobal_ (false)
    , fillComplete_ (false)
    , indicesAreSorted_ (true)
    , noRedundancies_ (true)
    , haveLocalConstants_ (false)
    , haveGlobalConstants_ (false)
    , haveRowInfo_ (true)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayRCP;
    using Teuchos::as;
    using Teuchos::rcp;
    typedef GlobalOrdinal GO;
    typedef LocalOrdinal LO;

    staticAssertions();
    const char tfecfFuncName[] = "CrsGraph(Map,Map,Kokkos::LocalStaticCrsGraph)";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! colMap.is_null (), std::runtime_error,
      ": The input column Map must be nonnull.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      k_local_graph_.numRows () != rowMap->getNodeNumElements (),
      std::runtime_error,
      ": The input row Map and the input local graph need to have the same "
      "number of rows.  The row Map claims " << rowMap->getNodeNumElements ()
      << " row(s), but the local graph claims " << k_local_graph_.numRows ()
      << " row(s).");
    // NOTE (mfh 17 Mar 2014) getNodeNumRows() returns
    // rowMap_->getNodeNumElements(), but it doesn't have to.
    // TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
    //   k_local_graph_.numRows () != getNodeNumRows (), std::runtime_error,
    //   ": The input row Map and the input local graph need to have the same "
    //   "number of rows.  The row Map claims " << getNodeNumRows () << " row(s), "
    //   "but the local graph claims " << k_local_graph_.numRows () << " row(s).");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! lclInds1D_.is_null () || ! gblInds1D_.is_null (), std::logic_error,
      ": cannot have 1D data structures allocated.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! lclInds2D_.is_null () || ! gblInds2D_.is_null (), std::logic_error,
      ": cannot have 2D data structures allocated.");

    nodeNumAllocated_ = k_local_graph_.row_map (getNodeNumRows ());
    nodeNumEntries_ = k_local_graph_.row_map (getNodeNumRows ());

    // NOTE (mfh 17 Mar 2014) We also need a version of this CrsGraph
    // constructor that takes a domain and range Map, as well as a row
    // and column Map.  In that case, we must pass the domain and
    // range Map into the following method.
    setDomainRangeMaps (rowMap_, rowMap_);
    makeImportExport ();

    RCP<ParameterList> lclparams;
    if (params == null) {
      lclparams = parameterList ();
    }
    else {
      lclparams = sublist (params, "Local Graph");
    }
    lclGraph_ = rcp (new local_graph_type (getRowMap ()->getNodeNumElements (),
                                           getColMap ()->getNodeNumElements (),
                                           getRowMap ()->getNode (), lclparams));
    ArrayRCP<const size_t> ptrs =
      arcp (k_lclGraph_.row_map.ptr_on_device (), 0,
            k_lclGraph_.row_map.dimension_0(),
            Kokkos::Compat::deallocator (k_lclGraph_.row_map), false);
    ArrayRCP<const LO> inds =
      arcp (k_lclGraph_.entries.ptr_on_device (), 0,
            k_lclGraph_.entries.dimension_0 (),
            Kokkos::Compat::deallocator (k_lclGraph_.entries), false);

    lclGraph_->setStructure (ptrs, inds);
    ptrs = null;
    inds = null;

    typename LocalStaticCrsGraphType::row_map_type d_ptrs = k_lclGraph_.row_map;
    typename LocalStaticCrsGraphType::entries_type d_inds = k_lclGraph_.entries;

    // Reset local properties
    upperTriangular_ = true;
    lowerTriangular_ = true;
    nodeMaxNumRowEntries_ = 0;
    nodeNumDiags_         = 0;

    // Compute triangular properties
    const size_t numLocalRows = getNodeNumRows ();
    for (size_t localRow = 0; localRow < numLocalRows; ++localRow) {
      const GO globalRow = rowMap_->getGlobalElement (localRow);
      const LO rlcid = colMap_->getLocalElement (globalRow);

      // It's entirely possible that the local matrix has no entries
      // in the column corresponding to the current row.  In that
      // case, the column Map may not necessarily contain that GID.
      // This is why we check whether rlcid is "invalid" (which means
      // that globalRow is not a GID in the column Map).
      if (rlcid != Teuchos::OrdinalTraits<LO>::invalid ()) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          rlcid >= as<LO> (d_ptrs.dimension_0 ()), std::runtime_error,
          ": The given row Map and/or column Map is/are not compatible with "
          "the provided local graph.");
        if (d_ptrs(rlcid) != d_ptrs(rlcid + 1)) {
          const size_t smallestCol = as<size_t> (d_inds(d_ptrs(rlcid)));
          const size_t largestCol = as<size_t> (d_inds(d_ptrs(rlcid + 1)));
          if (smallestCol < localRow) {
            upperTriangular_ = false;
          }
          if (localRow < largestCol) {
            lowerTriangular_ = false;
          }
          for (size_t i = d_ptrs(rlcid); i < d_ptrs(rlcid + 1); ++i) {
            if (d_inds(i) == rlcid) {
              ++nodeNumDiags_;
            }
          }
        }
        nodeMaxNumRowEntries_ =
          std::max (d_ptrs(rlcid + 1) - d_ptrs(rlcid), nodeMaxNumRowEntries_);
      }
    }

    haveLocalConstants_ = true;
    computeGlobalConstants();

    // finalize local graph
    Teuchos::EDiag diag = ( getNodeNumDiags() < getNodeNumRows() ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG );
    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    if      (isUpperTriangular()) uplo = Teuchos::UPPER_TRI;
    else if (isLowerTriangular()) uplo = Teuchos::LOWER_TRI;
    LocalMatOps::finalizeGraph(uplo,diag,*lclGraph_,params);
    fillComplete_ = true;

    checkInternalState();
  }
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::~CrsGraph()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  getValidParameters () const
  {
    RCP<ParameterList> params = parameterList ("Tpetra::CrsGraph");

    // Make a sublist for the Import.
    RCP<ParameterList> importSublist = parameterList ("Import");

    // FIXME (mfh 02 Apr 2012) We should really have the Import and
    // Export objects fill in these lists.  However, we don't want to
    // create an Import or Export unless we need them.  For now, we
    // know that the Import and Export just pass the list directly to
    // their Distributor, so we can create a Distributor here
    // (Distributor's constructor is a lightweight operation) and have
    // it fill in the list.

    // Fill in Distributor default parameters by creating a
    // Distributor and asking it to do the work.
    Distributor distributor (rowMap_->getComm(), importSublist);
    params->set ("Import", *importSublist, "How the Import performs communication.");

    // Make a sublist for the Export.  For now, it's a clone of the
    // Import sublist.  It's not a shallow copy, though, since we
    // might like the Import to do communication differently than the
    // Export.
    params->set ("Export", *importSublist, "How the Export performs communication.");

    return params;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  setParameterList (const RCP<ParameterList>& params)
  {
    RCP<const ParameterList> validParams = getValidParameters ();
    params->validateParametersAndSetDefaults (*validParams);
    this->setMyParamList (params);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalNumRows() const
  {
    return rowMap_->getGlobalNumElements();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalNumCols() const
  {
    const char tfecfFuncName[] = "getGlobalNumCols()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete() == false, std::runtime_error,
      ": requires domain map, which requires that fillComplete() has been "
      "called.");
    return getDomainMap()->getGlobalNumElements();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeNumRows() const
  {
    return rowMap_->getNodeNumElements();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeNumCols() const
  {
    const char tfecfFuncName[] = "getNodeNumCols()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasColMap() == false, std::runtime_error,
      ": requires column map.  This requires either that a custom column Map "
      "was given to the constructor, or that fillComplete() has been called.");
    return colMap_->getNodeNumElements();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeNumDiags() const
  {
    return nodeNumDiags_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalNumDiags() const
  {
    return globalNumDiags_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNode() const
  {
    return rowMap_->getNode();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getRowMap() const
  {
    return rowMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getColMap() const
  {
    return colMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getDomainMap() const
  {
    return domainMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getRangeMap() const
  {
    return rangeMap_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getImporter() const
  {
    return importer_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Export<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getExporter() const
  {
    return exporter_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::hasColMap() const
  {
    return (colMap_ != null);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isStorageOptimized() const
  {
    bool isOpt = indicesAreAllocated_ && (numRowEntries_ == null) && (getNodeNumRows() > 0);
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "isStorageOptimized()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (isOpt == true) && (getProfileType() == DynamicProfile), std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return isOpt;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ProfileType CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getProfileType() const
  {
    return pftype_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalNumEntries() const
  {
    return globalNumEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeNumEntries() const
  {
    return nodeNumEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  global_size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalMaxNumRowEntries() const
  {
    return globalMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeMaxNumRowEntries() const
  {
    return nodeMaxNumRowEntries_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isFillComplete() const
  {
    return fillComplete_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isFillActive() const
  {
    return !fillComplete_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isUpperTriangular() const
  {
    return upperTriangular_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isLowerTriangular() const
  {
    return lowerTriangular_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isLocallyIndexed() const
  {
    return indicesAreLocal_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isGloballyIndexed() const
  {
    return indicesAreGlobal_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeAllocationSize() const
  {
    return nodeNumAllocated_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const Comm<int> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getComm() const
  {
    return rowMap_->getComm();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  GlobalOrdinal CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getIndexBase() const
  {
    return rowMap_->getIndexBase();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::indicesAreAllocated() const
  {
    return indicesAreAllocated_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                    Internal utility methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isSorted() const
  {
    return indicesAreSorted_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::isMerged() const
  {
    return noRedundancies_;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::setLocallyModified ()
  {
    // FIXME (mfh 07 May 2013) How do we know that the change
    // introduced a redundancy, or even that it invalidated the sorted
    // order of indices?  CrsGraph has always made this conservative
    // guess.  It could be a bit costly to check at insertion time,
    // though.
    indicesAreSorted_ = false;
    noRedundancies_ = false;

    // We've modified the graph, so we'll have to recompute local
    // constants like the number of diagonal entries on this process.
    haveLocalConstants_ = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  allocateIndices (ELocalGlobal lg)
  {
    // This is a protected function, only callable by us.  If it was
    // called incorrectly, it is our fault.  That's why the tests
    // below throw std::logic_error instead of std::invalid_argument.
    const char tfecfFuncName[] = "allocateIndices()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed() && lg==GlobalIndices, std::logic_error,
      ": The graph is locally indexed, but Tpetra code is calling this method "
      "with lg=GlobalIndices.  Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed() && lg==LocalIndices, std::logic_error,
      ": The graph is globally indexed, but Tpetra code is calling this method "
      "with lg=LocalIndices.  Please report this bug to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated(), std::logic_error, ": The graph's indices are "
      "already allocated, but Tpetra code is calling allocateIndices() again.  "
      "Please report this bug to the Tpetra developers.");

    const size_t numRows = getNodeNumRows();
    indicesAreLocal_  = (lg == LocalIndices);
    indicesAreGlobal_ = (lg == GlobalIndices);
    nodeNumAllocated_ = 0;
    if (numAllocPerRow_ == null && getNodeNumRows() > 0) {
      // this wastes memory, temporarily, but it simplifies the code
      // and interfaces to follow
      ArrayRCP<size_t> tmpnumallocperrow = arcp<size_t>(numRows);
      std::fill(tmpnumallocperrow.begin(), tmpnumallocperrow.end(), numAllocForAllRows_);
      numAllocPerRow_ = tmpnumallocperrow;
    }
    //
    if (getProfileType() == StaticProfile) {
      //
      //  STATIC ALLOCATION PROFILE
      //
      // Have the local sparse kernels object allocate row offsets for
      // us, with first-touch allocation if applicable.  This is not
      // as important for global indices, because we never use global
      // indices in sparse kernels, but we might as well use the code
      // that we have for both the local and global indices cases.
      // Furthermore, first-touch allocation ensures that we don't
      // take up too much memory in any one NUMA affinity region.
      // rowPtrs_ = LocalMatOps::allocRowPtrs( getRowMap()->getNode(), numAllocPerRow_() );
      //std::cout<<"NumAllocPerRow: "<< numAllocPerRow_.size() << " " << numAllocPerRow_[0] << std::endl;
      k_rowPtrs_ = t_RowPtrs("Tpetra::CrsGraph::RowPtrs",numAllocPerRow_.size()+1);
      rowPtrs_ = Teuchos::arcp(k_rowPtrs_.ptr_on_device(), 0, k_rowPtrs_.dimension_0(),
                                         Kokkos::Compat::deallocator(k_rowPtrs_), false);

      // hack until we get parallel_scan in kokkos
      for(int i = 0; i < numAllocPerRow_.size(); i++) {
        k_rowPtrs_(i+1) = k_rowPtrs_(i)+numAllocPerRow_[i];
      }

      if (lg == LocalIndices) {
        // lclInds1D_ = LocalMatOps::template allocStorage<LocalOrdinal>( getRowMap()->getNode(), rowPtrs_() );
        k_lclInds1D_ = t_LocalOrdinal_1D("Tpetra::CrsGraph::lclInds1D_",*(rowPtrs_.end()-1));
        lclInds1D_ = Teuchos::arcp(k_lclInds1D_.ptr_on_device(), 0, k_lclInds1D_.dimension_0(),
                                   Kokkos::Compat::deallocator(k_lclInds1D_), false);
      }
      else {
        // gblInds1D_ = LocalMatOps::template allocStorage<GlobalOrdinal>( getRowMap()->getNode(), rowPtrs_() );
        k_gblInds1D_ = t_GlobalOrdinal_1D("Tpetra::CrsGraph::gblInds1D_",*(rowPtrs_.end()-1));
        gblInds1D_ = Teuchos::arcp(k_gblInds1D_.ptr_on_device(), 0, k_gblInds1D_.dimension_0(),
                                   Kokkos::Compat::deallocator(k_gblInds1D_), false);
      }
      nodeNumAllocated_ = rowPtrs_[numRows];
    }
    else {
      //
      //  DYNAMIC ALLOCATION PROFILE
      //
      typename ArrayRCP<const size_t>::iterator numalloc = numAllocPerRow_.begin();
      size_t howmany = numAllocForAllRows_;
      if (lg == LocalIndices) {
        lclInds2D_ = arcp< Array<LocalOrdinal> >(numRows);
        nodeNumAllocated_ = 0;
        for (size_t i=0; i < numRows; ++i) {
          if (numAllocPerRow_ != null) howmany = *numalloc++;
          nodeNumAllocated_ += howmany;
          if (howmany > 0) lclInds2D_[i].resize(howmany);
        }
      }
      else { // allocate global indices
        gblInds2D_ = arcp< Array<GlobalOrdinal> >(numRows);
        nodeNumAllocated_ = 0;
        for (size_t i=0; i < numRows; ++i) {
          if (numAllocPerRow_ != null) howmany = *numalloc++;
          nodeNumAllocated_ += howmany;
          if (howmany > 0) gblInds2D_[i].resize(howmany);
        }
      }
    }
    if (numRows > 0) {
      k_numRowEntries_ = t_numRowEntries_("Tpetra::CrsGraph::numRowEntries",numRows);
      numRowEntries_ = Teuchos::arcp(k_numRowEntries_.h_view.ptr_on_device(), 0, numRows,
                                     Kokkos::Compat::deallocator(k_numRowEntries_.h_view), false);
      //numRowEntries_ = arcp<size_t>(numRows);
      std::fill( numRowEntries_.begin(), numRowEntries_.end(), 0 );
    }
    // done with these
    numAllocForAllRows_ = 0;
    numAllocPerRow_     = null;
    indicesAreAllocated_ = true;
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  ArrayRCP<T> CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::allocateValues1D() const
  {
    const char tfecfFuncName[] = "allocateValues()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        indicesAreAllocated() == false,
        std::runtime_error, ": graph indices must be allocated before values."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() != StaticProfile,
        std::runtime_error, ": graph indices must be allocated in a static profile."
    );
    ArrayRCP<T> values1D;
    values1D = LocalMatOps::template allocStorage<T>( getRowMap()->getNode(), rowPtrs_() );
    return values1D;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class T>
  ArrayRCP<Array<T> >
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::allocateValues2D() const
  {
    const char tfecfFuncName[] = "allocateValues2D()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        indicesAreAllocated() == false,
        std::runtime_error, ": graph indices must be allocated before values."
    );
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() != DynamicProfile,
        std::runtime_error, ": graph indices must be allocated in a dynamic profile."
    );
    ArrayRCP<Array<T> > values2D;
    values2D = arcp<Array<T> > (getNodeNumRows ());
    if (lclInds2D_ != null) {
      const size_t numRows = lclInds2D_.size ();
      for (size_t r = 0; r < numRows; ++r) {
        values2D[r].resize (lclInds2D_[r].size ());
      }
    }
    else if (gblInds2D_ != null) {
      const size_t numRows = gblInds2D_.size ();
      for (size_t r = 0; r < numRows; ++r) {
        values2D[r].resize (gblInds2D_[r].size ());
      }
    }
    return values2D;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayView<const LocalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  getLocalView (RowInfo rowinfo) const
  {
    ArrayView<const LocalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (lclInds1D_ != null) {
        view = lclInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (!lclInds2D_[rowinfo.localRow].empty()) {
        view = lclInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayView<LocalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  getLocalViewNonConst (RowInfo rowinfo)
  {
    ArrayView<LocalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (lclInds1D_ != null) {
        view = lclInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (! lclInds2D_[rowinfo.localRow].empty()) {
        view = lclInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayView<const GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  getGlobalView (RowInfo rowinfo) const
  {
    ArrayView<const GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (! gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayView<GlobalOrdinal>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  getGlobalViewNonConst (RowInfo rowinfo)
  {
    ArrayView<GlobalOrdinal> view;
    if (rowinfo.allocSize > 0) {
      if (gblInds1D_ != null) {
        view = gblInds1D_ (rowinfo.offset1D, rowinfo.allocSize);
      }
      else if (!gblInds2D_[rowinfo.localRow].empty()) {
        view = gblInds2D_[rowinfo.localRow] ();
      }
    }
    return view;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RowInfo CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getRowInfo(size_t myRow) const
  {
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "getRowInfo()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        rowMap_->isNodeLocalElement (myRow) == false,
        std::logic_error,
        ": The given (local) row index myRow = " << myRow
        << " does not belong to the graph's row Map.  "
        "This probably indicates a bug in Tpetra::CrsGraph or Tpetra::CrsMatrix.  "
        "Please report this to the Tpetra developers.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasRowInfo() == false, std::logic_error,
      ": Late catch! Graph does not have row info anymore.  "
      "Error should have been caught earlier.  Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
    const size_t STINV = Teuchos::OrdinalTraits<size_t>::invalid ();
    RowInfo ret;
    ret.localRow = myRow;
    if (nodeNumAllocated_ != 0 && nodeNumAllocated_ != STINV) {
      // graph data structures have the info that we need
      //
      // if static graph, offsets tell us the allocation size
      if (getProfileType() == StaticProfile) {
        ret.offset1D  = rowPtrs_[myRow];
        ret.allocSize = rowPtrs_[myRow+1] - rowPtrs_[myRow];
        if (numRowEntries_ == null) {
          ret.numEntries = ret.allocSize;
        }
        else {
          ret.numEntries = numRowEntries_[myRow];
        }
      }
      else {
        ret.offset1D = STINV;
        if (isLocallyIndexed ()) {
          ret.allocSize = lclInds2D_[myRow].size();
        }
        else {
          ret.allocSize = gblInds2D_[myRow].size();
        }
        ret.numEntries = numRowEntries_[myRow];
      }
    }
    else if (nodeNumAllocated_ == 0) {
      // have performed allocation, but the graph has no allocation or entries
      ret.allocSize = 0;
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else if (indicesAreAllocated () == false) {
      // haven't performed allocation yet; probably won't hit this code
      if (numAllocPerRow_ == null) {
        ret.allocSize = numAllocForAllRows_;
      }
      else {
        ret.allocSize = numAllocPerRow_[myRow];
      }
      ret.numEntries = 0;
      ret.offset1D = STINV;
    }
    else {
      // don't know how we ended up here...
      TEUCHOS_TEST_FOR_EXCEPT(true);
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::staticAssertions() const
  {
    // Assumption: sizeof(GlobalOrdinal) >= sizeof(LocalOrdinal):
    //     This is so that we can store local indices in the memory
    //     formerly occupied by global indices.
    //
    // Assumption: max(GlobalOrdinal) >= max(LocalOrdinal) and
    //   max(size_t) >= max(LocalOrdinal)
    //     This is so that we can represent any LocalOrdinal as a
    //     size_t, and any LocalOrdinal as a GlobalOrdinal
    Teuchos::CompileTimeAssert<sizeof(GlobalOrdinal) < sizeof(LocalOrdinal)> cta_size1; (void)cta_size1;
    Teuchos::CompileTimeAssert<sizeof(global_size_t) < sizeof(size_t)      > cta_size2; (void)cta_size2;
    // can't call max() with CompileTimeAssert, because it isn't a constant expression; will need to make this a runtime check
    std::string msg = typeName(*this) + ": Object cannot be allocated with stated template arguments: size assumptions are not valid.";
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<LocalOrdinal>::max() > OrdinalTraits<size_t>::max(),          std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( (global_size_t)OrdinalTraits<LocalOrdinal>::max() > (global_size_t)OrdinalTraits<GlobalOrdinal>::max(),           std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( (size_t)OrdinalTraits<GlobalOrdinal>::max() > OrdinalTraits<global_size_t>::max(),  std::runtime_error, msg);
    TEUCHOS_TEST_FOR_EXCEPTION( OrdinalTraits<size_t>::max() > OrdinalTraits<global_size_t>::max(),                 std::runtime_error, msg);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertIndices (const RowInfo& rowinfo,
                 const SLocalGlobalViews &newInds,
                 const ELocalGlobal lg,
                 const ELocalGlobal I)
  {
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(
      lg != GlobalIndices && lg != LocalIndices, std::invalid_argument,
      "Tpetra::CrsGraph::insertIndices: lg must be either GlobalIndices or "
      "LocalIndices.");
#endif // HAVE_TPETRA_DEBUG
    size_t numNewInds = 0;
    if (lg == GlobalIndices) { // input indices are global
      ArrayView<const GlobalOrdinal> new_ginds = newInds.ginds;
      numNewInds = new_ginds.size();
      if (I == GlobalIndices) { // store global indices
        ArrayView<GlobalOrdinal> gind_view = getGlobalViewNonConst(rowinfo);
        std::copy(new_ginds.begin(), new_ginds.end(), gind_view.begin()+rowinfo.numEntries);
      }
      else if (I == LocalIndices) { // store local indices
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        typename ArrayView<const GlobalOrdinal>::iterator         in = new_ginds.begin();
        const typename ArrayView<const GlobalOrdinal>::iterator stop = new_ginds.end();
        typename ArrayView<LocalOrdinal>::iterator out = lind_view.begin()+rowinfo.numEntries;
        while (in != stop) {
          *out++ = colMap_->getLocalElement (*in++);
        }
      }
    }
    else if (lg == LocalIndices) { // input indices are local
      ArrayView<const LocalOrdinal> new_linds = newInds.linds;
      numNewInds = new_linds.size();
      if (I == LocalIndices) { // store local indices
        ArrayView<LocalOrdinal> lind_view = getLocalViewNonConst(rowinfo);
        std::copy(new_linds.begin(), new_linds.end(), lind_view.begin()+rowinfo.numEntries);
      }
      else if (I == GlobalIndices) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, "Tpetra::CrsGraph::"
          "insertIndices: the case where the input indices are local and the "
          "indices to write are global (lg=LocalIndices, I=GlobalIndices) is "
          "not implemented, because it does not make sense." << std::endl <<
          "If you have correct local column indices, that means the graph has "
          "a column Map.  In that case, you should be storing local indices.");
      }
    }
    numRowEntries_[rowinfo.localRow] += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();
    return numNewInds;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertGlobalIndicesImpl (const LocalOrdinal myRow,
                           const ArrayView<const GlobalOrdinal> &indices)
  {
    const char tfecfFuncName[] = "insertGlobalIndicesImpl";

    RowInfo rowInfo = getRowInfo(myRow);
    const size_t numNewInds = indices.size();
    const size_t newNumEntries = rowInfo.numEntries + numNewInds;
    if (newNumEntries > rowInfo.allocSize) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() == StaticProfile, std::runtime_error,
        ": new indices exceed statically allocated graph structure.");

      // update allocation, doubling size to reduce # reallocations
      size_t newAllocSize = 2*rowInfo.allocSize;
      if (newAllocSize < newNumEntries) {
        newAllocSize = newNumEntries;
      }
      gblInds2D_[myRow].resize(newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);
    }

    // Copy new indices at end of global index array
    if (gblInds1D_ != null)
      std::copy(indices.begin(), indices.end(),
                gblInds1D_.begin()+rowInfo.offset1D+rowInfo.numEntries);
    else
      std::copy(indices.begin(), indices.end(),
                gblInds2D_[myRow].begin()+rowInfo.numEntries);
    numRowEntries_[myRow] += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();

#ifdef HAVE_TPETRA_DEBUG
    {
      const size_t chkNewNumEntries = getNumEntriesInLocalRow (myRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        chkNewNumEntries != newNumEntries, std::logic_error,
        ": Internal logic error. Please contact Tpetra team.");
    }
#endif
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertLocalIndicesImpl (const LocalOrdinal myRow,
                          const ArrayView<const LocalOrdinal> &indices)
  {
    const char* tfecfFuncName ("insertLocallIndicesImpl");

    RowInfo rowInfo = getRowInfo(myRow);
    const size_t numNewInds = indices.size();
    const size_t newNumEntries = rowInfo.numEntries + numNewInds;
    if (newNumEntries > rowInfo.allocSize) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        getProfileType() == StaticProfile, std::runtime_error,
        ": new indices exceed statically allocated graph structure.");

      // update allocation, doubling size to reduce number of reallocations
      size_t newAllocSize = 2*rowInfo.allocSize;
      if (newAllocSize < newNumEntries)
        newAllocSize = newNumEntries;
      lclInds2D_[myRow].resize(newAllocSize);
      nodeNumAllocated_ += (newAllocSize - rowInfo.allocSize);
    }

    // Insert new indices at end of lclInds array
    if (lclInds1D_ != null)
      std::copy(indices.begin(), indices.end(),
                lclInds1D_.begin()+rowInfo.offset1D+rowInfo.numEntries);
    else
      std::copy(indices.begin(), indices.end(),
                lclInds2D_[myRow].begin()+rowInfo.numEntries);
    numRowEntries_[myRow] += numNewInds;
    nodeNumEntries_ += numNewInds;
    setLocallyModified ();
#ifdef HAVE_TPETRA_DEBUG
    {
      const size_t chkNewNumEntries = getNumEntriesInLocalRow (myRow);
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        chkNewNumEntries != newNumEntries, std::logic_error,
        ": Internal logic error. Please contact Tpetra team.");
    }
#endif
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class Scalar>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertIndicesAndValues (const RowInfo& rowInfo,
                          const SLocalGlobalViews& newInds,
                          const ArrayView<Scalar>& oldRowVals,
                          const ArrayView<const Scalar>& newRowVals,
                          const ELocalGlobal lg,
                          const ELocalGlobal I)
  {
    const size_t numNewInds = insertIndices (rowInfo, newInds, lg, I);
    typename ArrayView<const Scalar>::const_iterator newRowValsBegin =
      newRowVals.begin ();
    std::copy (newRowValsBegin, newRowValsBegin + numNewInds,
               oldRowVals.begin () + rowInfo.numEntries);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::sortRowIndices(RowInfo rowinfo)
  {
    if (rowinfo.numEntries > 0) {
      ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
      std::sort(inds_view.begin(), inds_view.begin() + rowinfo.numEntries);
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // in the future, perhaps this could use std::sort with a boost::zip_iterator
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template <class Scalar>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::sortRowIndicesAndValues(RowInfo rowinfo, ArrayView<Scalar> values)
  {
    if (rowinfo.numEntries > 0) {
      ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
      sort2(inds_view.begin(), inds_view.begin()+rowinfo.numEntries, values.begin());
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::mergeRowIndices(RowInfo rowinfo)
  {
    const char tfecfFuncName[] = "mergRowIndices()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized() == true, std::logic_error,
      ": The graph is already storage optimized, so we shouldn't be merging any indices."
      " Please report this bug to the Tpetra developers.");
    ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst(rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = std::unique(beg,end);
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the assignment below will destory the packed structure
    TEUCHOS_TEST_FOR_EXCEPT( isStorageOptimized() && mergedEntries != rowinfo.numEntries );
#endif
    numRowEntries_[rowinfo.localRow] = mergedEntries;
    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // in the future, this could use std::unique with a boost::zip_iterator
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  template<class Scalar>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  mergeRowIndicesAndValues (RowInfo rowinfo, const ArrayView<Scalar>& rowValues)
  {
    const char tfecfFuncName[] = "mergeRowIndicesAndValues";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized(), std::logic_error, ": It is invalid to call this "
      "method if the graph's storage has already been optimized." << std::endl
      << "Please report this bug to the Tpetra developers.");

    typedef typename ArrayView<Scalar>::iterator Iter;
    Iter rowValueIter = rowValues.begin ();
    ArrayView<LocalOrdinal> inds_view = getLocalViewNonConst (rowinfo);
    typename ArrayView<LocalOrdinal>::iterator beg, end, newend;

    // beg,end define a half-exclusive interval over which to iterate.
    beg = inds_view.begin();
    end = inds_view.begin() + rowinfo.numEntries;
    newend = beg;
    if (beg != end) {
      typename ArrayView<LocalOrdinal>::iterator cur = beg + 1;
      Iter vcur = rowValueIter + 1;
      Iter vend = rowValueIter;
      cur = beg+1;
      while (cur != end) {
        if (*cur != *newend) {
          // new entry; save it
          ++newend;
          ++vend;
          (*newend) = (*cur);
          (*vend) = (*vcur);
        }
        else {
          // old entry; merge it
          //(*vend) = f (*vend, *vcur);
          (*vend) += *vcur;
        }
        ++cur;
        ++vcur;
      }
      ++newend; // one past the last entry, per typical [beg,end) semantics
    }
    const size_t mergedEntries = newend - beg;
#ifdef HAVE_TPETRA_DEBUG
    // merge should not have eliminated any entries; if so, the
    // assignment below will destroy the packed structure
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isStorageOptimized() && mergedEntries != rowinfo.numEntries,
      std::logic_error,
      ": Merge was incorrect; it eliminated entries from the graph.  "
      << std::endl << "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    numRowEntries_[rowinfo.localRow] = mergedEntries;
    nodeNumEntries_ -= (rowinfo.numEntries - mergedEntries);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  setDomainRangeMaps (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& domainMap,
                      const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > >& rangeMap)
  {
    // simple pointer comparison for equality
    if (domainMap_ != domainMap) {
      domainMap_ = domainMap;
      importer_ = null;
    }
    if (rangeMap_ != rangeMap) {
      rangeMap_  = rangeMap;
      exporter_ = null;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  findLocalIndex (RowInfo rowinfo, LocalOrdinal ind, size_t hint) const
  {
    ArrayView<const LocalOrdinal> colInds = getLocalView (rowinfo);
    return this->findLocalIndex (rowinfo, ind, colInds, hint);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal,
            class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  findLocalIndex (RowInfo rowinfo,
                  LocalOrdinal ind,
                  ArrayView<const LocalOrdinal> colInds,
                  size_t hint) const
  {
    typedef typename ArrayView<const LocalOrdinal>::iterator IT;

    // If the hint was correct, then the hint is the offset to return.
    if (hint < rowinfo.numEntries && colInds[hint] == ind) {
      return hint;
    }

    // The hint was wrong, so we must search for the given column
    // index in the column indices for the given row.  How we do the
    // search depends on whether the graph's column indices are
    // sorted.
    IT beg = colInds.begin ();
    IT end = beg + rowinfo.numEntries;
    IT ptr = beg + rowinfo.numEntries; // "null"
    bool found = true;

    if (isSorted ()) {
      // binary search
      std::pair<IT,IT> p = std::equal_range (beg, end, ind);
      if (p.first == p.second) {
        found = false;
      }
      else {
        ptr = p.first;
      }
    }
    else {
      // direct search
      ptr = std::find (beg, end, ind);
      if (ptr == end) {
        found = false;
      }
    }

    if (found) {
      return Teuchos::as<size_t> (ptr - beg);
    }
    else {
      return Teuchos::OrdinalTraits<size_t>::invalid ();
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  findGlobalIndex (RowInfo rowinfo, GlobalOrdinal ind, size_t hint) const
  {
    typedef typename ArrayView<const GlobalOrdinal>::iterator IT;
    ArrayView<const GlobalOrdinal> indices = getGlobalView (rowinfo);

    // We don't actually require that the hint be a valid index.
    // If it is not in range, we just ignore it.
    if (hint < rowinfo.numEntries && indices[hint] == ind) {
      return hint;
    }

    IT beg = indices.begin ();
    IT end = indices.begin () + rowinfo.numEntries; // not indices.end()
    if (isSorted ()) { // use binary search
      const std::pair<IT,IT> p = std::equal_range (beg, end, ind);
      if (p.first == p.second) { // range of matching entries is empty
        return Teuchos::OrdinalTraits<size_t>::invalid ();
      } else {
        return p.first - beg;
      }
    }
    else { // not sorted; must use linear search
      const IT loc = std::find (beg, end, ind);
      if (loc == end) {
        return Teuchos::OrdinalTraits<size_t>::invalid ();
      } else {
        return loc - beg;
      }
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::clearGlobalConstants()
  {
    globalNumEntries_       = OrdinalTraits<global_size_t>::invalid();
    globalNumDiags_         = OrdinalTraits<global_size_t>::invalid();
    globalMaxNumRowEntries_ = OrdinalTraits<global_size_t>::invalid();
    haveGlobalConstants_    = false;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::checkInternalState() const
  {
#ifdef HAVE_TPETRA_DEBUG
    const global_size_t GSTI = OrdinalTraits<global_size_t>::invalid();
    const size_t         STI = OrdinalTraits<size_t>::invalid();
    std::string err = typeName(*this) + "::checkInternalState(): Likely internal logic error. Please contact Tpetra team.";
    // check the internal state of this data structure
    // this is called by numerous state-changing methods, in a debug build, to ensure that the object
    // always remains in a valid state
    // the graph should have been allocated with a row map
    TEUCHOS_TEST_FOR_EXCEPTION( rowMap_ == null,                     std::logic_error, err );
    // am either complete or active
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() == isFillComplete(),  std::logic_error, err );
    // if active, i have no local graph
    TEUCHOS_TEST_FOR_EXCEPTION( isFillActive() && lclGraph_ != null, std::logic_error, err );
    // if the graph has been fill completed, then all maps should be present
    TEUCHOS_TEST_FOR_EXCEPTION( isFillComplete() == true && (colMap_ == null || rangeMap_ == null || domainMap_ == null), std::logic_error, err );
    // if storage has been optimized, then indices should have been allocated (even if trivially so)
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && indicesAreAllocated() == false, std::logic_error, err );
    // if storage has been optimized, then number of allocated is now the number of entries
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() == true && getNodeAllocationSize() != getNodeNumEntries(), std::logic_error, err );
    // if graph doesn't have the global constants, then they should all be marked as invalid
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == false && ( globalNumEntries_ != GSTI || globalNumDiags_ != GSTI || globalMaxNumRowEntries_ != GSTI ), std::logic_error, err );
    // if the graph has global cosntants, then they should be valid.
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ == GSTI || globalNumDiags_ == GSTI || globalMaxNumRowEntries_ == GSTI ), std::logic_error, err );
    TEUCHOS_TEST_FOR_EXCEPTION( haveGlobalConstants_ == true && ( globalNumEntries_ < nodeNumEntries_ || globalNumDiags_ < nodeNumDiags_ || globalMaxNumRowEntries_ < nodeMaxNumRowEntries_ ),
                        std::logic_error, err );
    // if indices are allocated, then the allocation specifications should have been released
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == true  && (numAllocForAllRows_ != 0 || numAllocPerRow_ != null),                        std::logic_error, err );
    // if indices are not allocated, then information dictating allocation quantities should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && (nodeNumAllocated_ != STI || nodeNumEntries_ != 0),                           std::logic_error, err );
    // if storage is optimized, then profile should be static
    TEUCHOS_TEST_FOR_EXCEPTION( isStorageOptimized() && pftype_ != StaticProfile,                                                               std::logic_error, err );
    // if rowPtrs_ exists, it is required to have N+1 rows, and rowPtrs_[N] == gblInds1D_.size()/lclInds1D_.size()
    TEUCHOS_TEST_FOR_EXCEPTION( isGloballyIndexed() && rowPtrs_ != null && ((size_t)rowPtrs_.size() != getNodeNumRows()+1 || rowPtrs_[getNodeNumRows()] != (size_t)gblInds1D_.size()), std::logic_error, err );
    TEUCHOS_TEST_FOR_EXCEPTION(  isLocallyIndexed() && rowPtrs_ != null && ((size_t)rowPtrs_.size() != getNodeNumRows()+1 || rowPtrs_[getNodeNumRows()] != (size_t)lclInds1D_.size()), std::logic_error, err );
    // if profile is dynamic and we have allocated, then 2D allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && indicesAreAllocated() && getNodeNumRows() > 0 && lclInds2D_ == null && gblInds2D_ == null,
                                                                                                                                        std::logic_error, err );
    // if profile is dynamic, then numentries and 2D indices are needed and should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && indicesAreAllocated() && getNodeNumRows() > 0 && (numRowEntries_ == null || (lclInds2D_ == null && gblInds2D_ == null)),
                                                                                                                                        std::logic_error, err );
    // if profile is dynamic, then 1D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && (lclInds1D_ != null || gblInds1D_ != null),                                std::logic_error, err );
    // if profile is dynamic, then row offsets should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == DynamicProfile && rowPtrs_ != null,                                                          std::logic_error, err );
    // if profile is static and we have allocated non-trivially, then 1D allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == StaticProfile && indicesAreAllocated() && getNodeAllocationSize() > 0 && lclInds1D_ == null && gblInds1D_ == null,
                                                                                                                                        std::logic_error, err );
    // if profile is static, then 2D allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( pftype_ == StaticProfile && (lclInds2D_ != null || gblInds2D_ != null),                                 std::logic_error, err );
    // if indices are not allocated, then none of the buffers should be.
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && (rowPtrs_ != null || numRowEntries_ != null ||
                                                                 lclInds1D_ != null || lclInds2D_ != null ||
                                                                 gblInds1D_ != null || gblInds2D_ != null),                             std::logic_error, err );
    // indices may be local or global only if they are allocated (numAllocated is redundant; could simply be indicesAreLocal_ || indicesAreGlobal_)
    TEUCHOS_TEST_FOR_EXCEPTION( (indicesAreLocal_ == true || indicesAreGlobal_ == true) && indicesAreAllocated_ == false,               std::logic_error, err );
    // indices may be local or global, but not both
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && indicesAreGlobal_ == true,                                                  std::logic_error, err );
    // if indices are local, then global allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && (gblInds1D_ != null || gblInds2D_ != null),                                 std::logic_error, err );
    // if indices are global, then local allocations should not be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && (lclInds1D_ != null || lclInds2D_ != null),                                std::logic_error, err );
    // if indices are local, then local allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreLocal_ == true && getNodeAllocationSize() > 0 && lclInds1D_ == null && getNodeNumRows() > 0 && lclInds2D_ == null,
                                                                                                                          std::logic_error, err );
    // if indices are global, then global allocations should be present
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreGlobal_ == true && getNodeAllocationSize() > 0 && gblInds1D_ == null && getNodeNumRows() > 0 && gblInds2D_ == null,
                                                                                                                          std::logic_error, err );
    // if indices are allocated, then we should have recorded how many were allocated
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == true  && nodeNumAllocated_ == STI,                                       std::logic_error, err );
    // if indices are not allocated, then the allocation size should be marked invalid
    TEUCHOS_TEST_FOR_EXCEPTION( indicesAreAllocated() == false && nodeNumAllocated_ != STI,                                       std::logic_error, err );
    // check the actual allocations
    if (indicesAreAllocated()) {
      size_t actualNumAllocated = 0;
      if (pftype_ == DynamicProfile) {
        if (isGloballyIndexed() && gblInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += gblInds2D_[r].size();
          }
        }
        else if (isLocallyIndexed() && lclInds2D_ != null) {
          for (size_t r = 0; r < getNodeNumRows(); ++r) {
            actualNumAllocated += lclInds2D_[r].size();
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
      }
      else if (rowPtrs_ != null) { // pftype_ == StaticProfile
        actualNumAllocated = rowPtrs_[getNodeNumRows()];
        TEUCHOS_TEST_FOR_EXCEPTION(  isLocallyIndexed() == true && (size_t)lclInds1D_.size() != actualNumAllocated, std::logic_error, err );
        TEUCHOS_TEST_FOR_EXCEPTION( isGloballyIndexed() == true && (size_t)gblInds1D_.size() != actualNumAllocated, std::logic_error, err );
        TEUCHOS_TEST_FOR_EXCEPTION(actualNumAllocated != nodeNumAllocated_, std::logic_error, err );
      }
    }
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  //                                                                         //
  //                  User-visible class methods                             //
  //                                                                         //
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNumEntriesInGlobalRow(GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (hasRowInfo() && lrow != OrdinalTraits<LocalOrdinal>::invalid())
    {
      RowInfo rowinfo = getRowInfo(lrow);
      ret = rowinfo.numEntries;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNumEntriesInLocalRow(LocalOrdinal localRow) const
  {
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (hasRowInfo() && rowMap_->isNodeLocalElement(localRow)) {
      RowInfo rowinfo = getRowInfo(localRow);
      ret = rowinfo.numEntries;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNumAllocatedEntriesInGlobalRow(GlobalOrdinal globalRow) const
  {
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (hasRowInfo() && lrow != OrdinalTraits<LocalOrdinal>::invalid())
    {
      RowInfo rowinfo = getRowInfo(lrow);
      ret = rowinfo.allocSize;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  size_t CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNumAllocatedEntriesInLocalRow(LocalOrdinal localRow) const
  {
    size_t ret = OrdinalTraits<size_t>::invalid();
    if (hasRowInfo() && rowMap_->isNodeLocalElement(localRow)) {
      RowInfo rowinfo = getRowInfo(localRow);
      ret = rowinfo.allocSize;
    }
    return ret;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayRCP<const size_t> CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodeRowPtrs() const
  {
    return rowPtrs_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  ArrayRCP<const LocalOrdinal> CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getNodePackedIndices() const
  {
    return lclInds1D_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getLocalRowCopy(LocalOrdinal localRow, const ArrayView<LocalOrdinal> &indices, size_t &NumIndices) const
  {
    // can only do this if
    // * we have local indices: isLocallyIndexed()
    // or
    // * we are capable of producing them: isGloballyIndexed() && hasColMap()
    // short circuit if we aren't allocated
    const char tfecfFuncName[] = "getLocalRowCopy(localRow,...)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true && hasColMap() == false, std::runtime_error, ": local indices cannot be produced.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
        ": localRow (== " << localRow << ") is not valid on this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(hasRowInfo() == false, std::runtime_error, ": graph row information was deleted at fillComplete().");
    const RowInfo rowinfo = getRowInfo(localRow);
    NumIndices = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)indices.size() < NumIndices, std::runtime_error,
        ": specified storage (size==" << indices.size()
        << ") is not large enough to hold all entries for this row (NumIndices == " << NumIndices << ").");
    if (isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> lview = getLocalView(rowinfo);
      std::copy( lview.begin(), lview.begin() + NumIndices, indices.begin());
    }
    else if (isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> gview = getGlobalView(rowinfo);
      for (size_t j=0; j < NumIndices; ++j) {
        indices[j] = colMap_->getLocalElement(gview[j]);
      }
    }
    else {
#ifdef HAVE_TPETRA_DEBUG
      // should have fallen in one of the above
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( indicesAreAllocated() == true, std::logic_error, ": Internal logic error. Please contact Tpetra team.");
#endif
      NumIndices = 0;
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalRowCopy(GlobalOrdinal globalRow, const ArrayView<GlobalOrdinal> &indices, size_t &NumIndices) const
  {
    // we either currently store global indices, or we have a column map with which to transcribe our local indices for the user
    const LocalOrdinal lrow = rowMap_->getLocalElement(globalRow);
    const char tfecfFuncName[] = "getGlobalRowCopy(globalRow,...)";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(lrow == OrdinalTraits<LocalOrdinal>::invalid(), std::runtime_error,
        ": globalRow (== " << globalRow << ") does not belong to this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(hasRowInfo() == false, std::runtime_error, ": graph row information was deleted at fillComplete().");
    const RowInfo rowinfo = getRowInfo((size_t)lrow);
    NumIndices = rowinfo.numEntries;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC((size_t)indices.size() < NumIndices, std::runtime_error,
        ": specified storage (size==" << indices.size()
        << ") is not large enough to hold all entries for this row (NumIndices == " << NumIndices << ").");
    if (isLocallyIndexed()) {
      ArrayView<const LocalOrdinal> lview = getLocalView(rowinfo);

      for (size_t j=0; j < NumIndices; ++j) {
        indices[j] = colMap_->getGlobalElement(lview[j]);
      }
    }
    else if (isGloballyIndexed()) {
      ArrayView<const GlobalOrdinal> gview = getGlobalView(rowinfo);
      std::copy(gview.begin(), gview.begin() + NumIndices, indices.begin());
    }
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getLocalRowView(LocalOrdinal localRow, ArrayView<const LocalOrdinal> &indices) const
  {
    const char tfecfFuncName[] = "getLocalRowView()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isGloballyIndexed() == true, std::runtime_error, ": local indices cannot be provided.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(hasRowInfo() == false, std::runtime_error, ": graph row information was deleted at fillComplete().");
    indices = null;
    if (rowMap_->isNodeLocalElement(localRow) == true) {
      const RowInfo rowinfo = getRowInfo(localRow);
      if (rowinfo.numEntries > 0) {
        indices = getLocalView(rowinfo);
        indices = indices(0,rowinfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)indices.size() != getNumEntriesInLocalRow(localRow), std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    return;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getGlobalRowView(GlobalOrdinal globalRow, ArrayView<const GlobalOrdinal> &indices) const
  {
    using Teuchos::as;
    const char tfecfFuncName[] = "getGlobalRowView()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      ": The graph is locally indexed, so we cannot return a view with global "
      "column indices.  Use getGlobalRowCopy() instead.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! hasRowInfo (), std::runtime_error,
      ": graph row information was deleted at fillComplete().");

    // isNodeGlobalElement() requires a global to local lookup anyway,
    // and getLocalElement() returns invalid() if the element wasn't found.
    const LocalOrdinal localRow = rowMap_->getLocalElement (globalRow);
    indices = null;
    if (localRow != Teuchos::OrdinalTraits<LocalOrdinal>::invalid ()) {
      const RowInfo rowInfo = getRowInfo (as<size_t> (localRow));
      if (rowInfo.numEntries > 0) {
        indices = (getGlobalView (rowInfo)) (0, rowInfo.numEntries);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    using std::endl;
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<size_t> (indices.size ()) != getNumEntriesInGlobalRow (globalRow),
      std::logic_error,
      ": Violated stated post-conditions:"
      << "  indices.size () = " << indices.size () << endl
      << "  as<size_t> (indices.size ()) = " << as<size_t> (indices.size ())
      << endl << "  getNumEntriesInGlobalRow (globalRow = " << globalRow
      << ") = " << getNumEntriesInGlobalRow (globalRow) << endl
      << "Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertLocalIndices (const LocalOrdinal localRow,
                      const ArrayView<const LocalOrdinal> &indices)
  {
    const char tfecfFuncName[] = "insertLocalIndices";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed() == true, std::runtime_error,
      ": graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasColMap() == false, std::runtime_error,
      ": cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
      ": row does not belong to this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasRowInfo() == false, std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    if (! indicesAreAllocated ()) {
      allocateIndices (LocalIndices);
    }

#ifdef HAVE_TPETRA_DEBUG
    // In a debug build, if the graph has a column Map, test whether
    // any of the given column indices are not in the column Map.
    // Keep track of the invalid column indices so we can tell the
    // user about them.
    if (hasColMap ()) {
      using Teuchos::Array;
      using Teuchos::toString;
      using std::endl;
      typedef typename Teuchos::ArrayView<const LocalOrdinal>::size_type size_type;

      const map_type& colMap = * (getColMap ());
      Array<LocalOrdinal> badColInds;
      bool allInColMap = true;
      for (size_type k = 0; k < indices.size (); ++k) {
        if (! colMap.isNodeLocalElement (indices[k])) {
          allInColMap = false;
          badColInds.push_back (indices[k]);
        }
      }
      if (! allInColMap) {
        std::ostringstream os;
        os << "Tpetra::CrsMatrix::insertLocalIndices: You attempted to insert "
          "entries in owned row " << localRow << ", at the following column "
          "indices: " << toString (indices) << "." << endl;
        os << "Of those, the following indices are not in the column Map on "
          "this process: " << toString (badColInds) << "." << endl << "Since "
          "the graph has a column Map already, it is invalid to insert entries "
          "at those locations.";
        TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
      }
    }
#endif // HAVE_TPETRA_DEBUG

    insertLocalIndicesImpl (localRow, indices);

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isLocallyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertLocalIndicesFiltered (const LocalOrdinal localRow,
                              const ArrayView<const LocalOrdinal> &indices)
  {
    typedef LocalOrdinal LO;
    const char tfecfFuncName[] = "insertLocalIndicesFiltered";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isGloballyIndexed() == true, std::runtime_error,
      ": graph indices are global; use insertGlobalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasColMap() == false, std::runtime_error,
      ": cannot insert local indices without a column map.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      rowMap_->isNodeLocalElement(localRow) == false, std::runtime_error,
      ": row does not belong to this node.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasRowInfo() == false, std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    if (indicesAreAllocated() == false) {
      allocateIndices(LocalIndices);
    }

     // If we have a column map, use it to filter the entries.
    if (hasColMap ()) {
      Array<LO> filtered_indices(indices);
      SLocalGlobalViews inds_view;
      SLocalGlobalNCViews inds_ncview;
      inds_ncview.linds = filtered_indices();
      const size_t numFilteredEntries =
        filterIndices<LocalIndices>(inds_ncview);
      inds_view.linds = filtered_indices (0, numFilteredEntries);
      insertLocalIndicesImpl(localRow, inds_view.linds);
    }
    else {
      insertLocalIndicesImpl(localRow, indices);
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isLocallyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertGlobalIndices (const GlobalOrdinal grow,
                       const ArrayView<const GlobalOrdinal> &indices)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;
    const char tfecfFuncName[] = "insertGlobalIndices";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed() == true, std::runtime_error,
      ": graph indices are local; use insertLocalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasRowInfo() == false, std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": You are not allowed to call this method if fill is not active.  "
      "If fillComplete has been called, you must first call resumeFill "
      "before you may insert indices.");
    if (! indicesAreAllocated ()) {
      allocateIndices (GlobalIndices);
    }
    const LO myRow = rowMap_->getLocalElement (grow);
    if (myRow != Teuchos::OrdinalTraits<LO>::invalid ()) {
#ifdef HAVE_TPETRA_DEBUG
      if (hasColMap ()) {
        using std::endl;
        const map_type& colMap = * (getColMap ());
        // In a debug build, keep track of the nonowned ("bad") column
        // indices, so that we can display them in the exception
        // message.  In a release build, just ditch the loop early if
        // we encounter a nonowned column index.
        Array<GO> badColInds;
        bool allInColMap = true;
        for (size_type k = 0; k < indices.size (); ++k) {
          if (! colMap.isNodeGlobalElement (indices[k])) {
            allInColMap = false;
            badColInds.push_back (indices[k]);
          }
        }
        if (! allInColMap) {
          std::ostringstream os;
          os << "Tpetra::CrsGraph::insertGlobalIndices: You attempted to insert "
            "entries in owned row " << grow << ", at the following column "
            "indices: " << toString (indices) << "." << endl;
          os << "Of those, the following indices are not in the column Map on "
            "this process: " << toString (badColInds) << "." << endl << "Since "
            "the matrix has a column Map already, it is invalid to insert "
            "entries at those locations.";
          TEUCHOS_TEST_FOR_EXCEPTION(! allInColMap, std::invalid_argument, os.str ());
        }
      }
#endif // HAVE_TPETRA_DEBUG
      insertGlobalIndicesImpl (myRow, indices);
    }
    else { // a nonlocal row
      const size_type numIndices = indices.size ();
      // This creates the Array if it doesn't exist yet.  std::map's
      // operator[] does a lookup each time, so it's better to pull
      // nonlocals_[grow] out of the loop.
      std::vector<GO>& nonlocalRow = nonlocals_[grow];
      for (size_type k = 0; k < numIndices; ++k) {
        nonlocalRow.push_back (indices[k]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isGloballyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  insertGlobalIndicesFiltered (const GlobalOrdinal grow,
                               const ArrayView<const GlobalOrdinal> &indices)
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "insertGlobalIndicesFiltered";

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed() == true, std::runtime_error,
      ": graph indices are local; use insertLocalIndices().");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      hasRowInfo() == false, std::runtime_error,
      ": graph row information was deleted at fillComplete().");
    // This can't really be satisfied for now, because if we are
    // fillComplete(), then we are local.  In the future, this may
    // change.  However, the rule that modification require active
    // fill will not change.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillActive() == false, std::runtime_error,
      ": You are not allowed to call this method if fill is not active.  "
      "If fillComplete has been called, you must first call resumeFill "
      "before you may insert indices.");
    if (indicesAreAllocated() == false) {
      allocateIndices (GlobalIndices);
    }
    const LO myRow = rowMap_->getLocalElement (grow);
    if (myRow != Teuchos::OrdinalTraits<LO>::invalid ()) {
      // If we have a column map, use it to filter the entries.
      if (hasColMap ()) {
        Array<GO> filtered_indices(indices);
        SLocalGlobalViews inds_view;
        SLocalGlobalNCViews inds_ncview;
        inds_ncview.ginds = filtered_indices();
        const size_t numFilteredEntries =
          filterIndices<GlobalIndices> (inds_ncview);
        inds_view.ginds = filtered_indices (0, numFilteredEntries);
        insertGlobalIndicesImpl(myRow, inds_view.ginds);
      }
      else {
       insertGlobalIndicesImpl(myRow, indices);
      }
    }
    else { // a nonlocal row
      typedef typename ArrayView<const GO>::size_type size_type;
      const size_type numIndices = indices.size ();
      for (size_type k = 0; k < numIndices; ++k) {
        nonlocals_[grow].push_back (indices[k]);
      }
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      indicesAreAllocated() == false || isGloballyIndexed() == false,
      std::logic_error,
      ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::removeLocalIndices(LocalOrdinal lrow)
  {
    const char tfecfFuncName[] = "removeLocalIndices()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == false,                    std::runtime_error, ": requires that fill is active.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isStorageOptimized() == true,               std::runtime_error, ": cannot remove indices after optimizeStorage() has been called.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isGloballyIndexed() == true,                std::runtime_error, ": graph indices are global.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( rowMap_->isNodeLocalElement(lrow) == false, std::runtime_error, ": row does not belong to this node.");
    if (indicesAreAllocated() == false) {
      allocateIndices(LocalIndices);
    }
    //
    clearGlobalConstants();
    //
    if (numRowEntries_ != null) {
      nodeNumEntries_ -= numRowEntries_[lrow];
      numRowEntries_[lrow] = 0;
    }
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(getNumEntriesInLocalRow(lrow) != 0 || indicesAreAllocated() == false || isLocallyIndexed() == false, std::logic_error,
        ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::setAllIndices(const t_RowPtrs & rowPointers,const t_LocalOrdinal_1D & columnIndices)
  {
    const char tfecfFuncName[] = "setAllIndices()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( hasColMap() == false, std::runtime_error, ": requires a ColMap.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)rowPointers.size() != getNodeNumRows()+1, std::runtime_error, ": requires rowPointers.size() == getNodeNumRows()+1");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( lclInds1D_ != Teuchos::null || gblInds1D_ != Teuchos::null, std::runtime_error, ": cannot have 1D data structures allocated.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( lclInds2D_ != Teuchos::null || gblInds2D_ != Teuchos::null, std::runtime_error, ": cannot have 2D data structures allocated.");

    indicesAreAllocated_ = true;
    indicesAreLocal_     = true;
    pftype_              = StaticProfile; // if the profile wasn't static before, it sure is now.
    k_lclInds1D_         = columnIndices;
    lclInds1D_           = Kokkos::Compat::persistingView(k_lclInds1D_);
    k_rowPtrs_           = rowPointers;
    rowPtrs_             = Kokkos::Compat::persistingView(k_rowPtrs_);
    nodeNumEntries_ = nodeNumAllocated_ = rowPtrs_[getNodeNumRows()];
    checkInternalState();
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::setAllIndices(const ArrayRCP<size_t> & rowPointers,const ArrayRCP<LocalOrdinal> & columnIndices)
  {
    const char tfecfFuncName[] = "setAllIndices()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( hasColMap() == false, std::runtime_error, ": requires a ColMap.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( (size_t)rowPointers.size() != getNodeNumRows()+1, std::runtime_error, ": requires rowPointers.size() == getNodeNumRows()+1");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( lclInds1D_ != Teuchos::null || gblInds1D_ != Teuchos::null, std::runtime_error, ": cannot have 1D data structures allocated.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( lclInds2D_ != Teuchos::null || gblInds2D_ != Teuchos::null, std::runtime_error, ": cannot have 2D data structures allocated.");

    indicesAreAllocated_ = true;
    indicesAreLocal_     = true;
    pftype_              = StaticProfile; // if the profile wasn't static before, it sure is now.
    k_lclInds1D_         = Kokkos::Compat::getKokkosViewDeepCopy<DeviceType>(columnIndices());
    lclInds1D_           = columnIndices;
    k_rowPtrs_           = Kokkos::Compat::getKokkosViewDeepCopy<DeviceType>(rowPointers());
    rowPtrs_             = rowPointers;
    nodeNumEntries_ = nodeNumAllocated_ = rowPtrs_[getNodeNumRows()];
    checkInternalState();
  }

  // TODO: in the future, globalAssemble() should use import/export functionality
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::globalAssemble()
  {
    using Teuchos::Array;
    using Teuchos::as;
    using Teuchos::Comm;
    using Teuchos::gatherAll;
    using Teuchos::ireceive;
    using Teuchos::isend;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    using Teuchos::toString;
    using Teuchos::waitAll;
    using std::endl;
    using std::make_pair;
    using std::pair;
    typedef GlobalOrdinal GO;
    typedef typename std::map<GO, std::vector<GO> >::const_iterator NLITER;
    typedef typename Array<GO>::size_type size_type;

    const char tfecfFuncName[] = "globalAssemble"; // for exception macro
    RCP<const Comm<int> > comm = getComm();

    const int numImages = comm->getSize();
    const int myImageID = comm->getRank();
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive(), std::runtime_error,
      ": requires that fill is active.");
    // Determine if any nodes have global entries to share
    {
      size_t MyNonlocals = nonlocals_.size(), MaxGlobalNonlocals;
      reduceAll (*comm, REDUCE_MAX, MyNonlocals, outArg (MaxGlobalNonlocals));
      if (MaxGlobalNonlocals == 0) {
        return;  // no entries to share
      }
    }

    // compute a list of NLRs from nonlocals_ and use it to compute:
    //      IdsAndRows: a vector of (id,row) pairs
    //          NLR2Id: a map from NLR to the Id that owns it
    // globalNeighbors: a global graph of connectivity between images:
    //   globalNeighbors(i,j) indicates that j sends to i
    //         sendIDs: a list of all images I send to
    //         recvIDs: a list of all images I receive from (constructed later)
    Array<pair<int, GO> > IdsAndRows;
    std::map<GO, int> NLR2Id;
    Teuchos::SerialDenseMatrix<int, char> globalNeighbors;
    Array<int> sendIDs, recvIDs;
    {
      // nonlocals_ contains the entries we are holding for all
      // nonowned rows.  Compute list of rows for which we have data.
      Array<GO> NLRs;
      std::set<GO> setOfRows;
      for (NLITER iter = nonlocals_.begin (); iter != nonlocals_.end (); ++iter) {
        setOfRows.insert (iter->first);
      }
      // copy the elements in the set into an Array
      NLRs.resize (setOfRows.size ());
      std::copy (setOfRows.begin (), setOfRows.end (), NLRs.begin ());

      // get a list of ImageIDs for the non-local rows (NLRs)
      Array<int> NLRIds(NLRs.size());
      {
        const LookupStatus stat =
          rowMap_->getRemoteIndexList (NLRs (), NLRIds ());
        int lclerror = ( stat == IDNotPresent ? 1 : 0 );
        int gblerror;
        reduceAll<int, int> (*comm, REDUCE_MAX, lclerror, outArg (gblerror));
        if (gblerror != 0) {
          const int myRank = comm->getRank ();
          std::ostringstream os;
          os << "On one or more processes in the communicator, "
             << "there were insertions into rows of the graph that do not "
             << "exist in the row Map on any process in the communicator."
             << endl << "This process " << myRank << " is "
             << (lclerror == 0 ? "not " : "") << "one of those offending "
             << "processes." << endl;
          if (lclerror != 0) {
            // If NLRIds[k] is -1, then NLRs[k] is a row index not in
            // the row Map.  Collect this list of invalid row indices
            // for display in the exception message.
            Array<GO> invalidNonlocalRows;
            for (size_type k = 0; k < NLRs.size (); ++k) {
              if (NLRIds[k] == -1) {
                invalidNonlocalRows.push_back (NLRs[k]);
              }
            }
            const size_type numInvalid = invalidNonlocalRows.size ();
            os << "On this process, " << numInvalid << " nonlocal row"
               << (numInvalid != 1 ? "s " : " ") << " were inserted that are "
               << "not in the row Map on any process." << endl;
            // Don't print _too_ many nonlocal rows.
            if (numInvalid <= 100) {
              os << "Offending row indices: "
                 << toString (invalidNonlocalRows ()) << endl;
            }
          }
          TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
            gblerror != 0, std::runtime_error,
            ": nonlocal entries correspond to invalid rows."
            << endl << os.str ());
        }
      }

      // build up a list of neighbors, as well as a map between NLRs and Ids
      // localNeighbors[i] != 0 iff I have data to send to image i
      // put NLRs,Ids into an array of pairs
      IdsAndRows.reserve(NLRs.size());
      Array<char> localNeighbors(numImages, 0);
      typename Array<GO>::const_iterator nlr;
      typename Array<int>::const_iterator id;
      for (nlr = NLRs.begin(), id = NLRIds.begin();
           nlr != NLRs.end(); ++nlr, ++id) {
        NLR2Id[*nlr] = *id;
        localNeighbors[*id] = 1;
        // IdsAndRows.push_back(make_pair<int,GlobalOrdinal>(*id,*nlr));
        IdsAndRows.push_back(make_pair(*id,*nlr));
      }
      for (int j=0; j<numImages; ++j) {
        if (localNeighbors[j]) {
          sendIDs.push_back(j);
        }
      }
      // sort IdsAndRows, by Ids first, then rows
      std::sort(IdsAndRows.begin(),IdsAndRows.end());
      // gather from other nodes to form the full graph
      globalNeighbors.shapeUninitialized(numImages,numImages);
      gatherAll (*getComm(), numImages, localNeighbors.getRawPtr(),
                 numImages * numImages, globalNeighbors.values());
      // globalNeighbors at this point contains (on all images) the
      // connectivity between the images.
      // globalNeighbors(i,j) != 0 means that j sends to i/that i receives from j
    }

    //////////////////////////////////////////////////////////////////////////////////////
    // FIGURE OUT WHO IS SENDING TO WHOM AND HOW MUCH
    // DO THIS IN THE PROCESS OF PACKING ALL OUTGOING DATA ACCORDING TO DESTINATION ID
    //////////////////////////////////////////////////////////////////////////////////////

    // loop over all columns to know from which images I can expect to receive something
    for (int j=0; j<numImages; ++j) {
      if (globalNeighbors(myImageID,j)) {
        recvIDs.push_back(j);
      }
    }
    const size_t numRecvs = recvIDs.size();

    // we know how many we're sending to already
    // form a contiguous list of all data to be sent
    // track the number of entries for each ID
    Array<pair<GO, GO> > IJSendBuffer;
    Array<size_t> sendSizes(sendIDs.size(), 0);
    size_t numSends = 0;
    for (typename Array<pair<int, GO> >::const_iterator IdAndRow = IdsAndRows.begin();
         IdAndRow != IdsAndRows.end(); ++IdAndRow) {
      int id = IdAndRow->first;
      GO row = IdAndRow->second;
      // have we advanced to a new send?
      if (sendIDs[numSends] != id) {
        numSends++;
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(sendIDs[numSends] != id,
          std::logic_error, ": internal logic error. Contact Tpetra team.");
      }
      // copy data for row into contiguous storage
      for (typename std::vector<GO>::const_iterator j = nonlocals_[row].begin ();
           j != nonlocals_[row].end (); ++j) {
        IJSendBuffer.push_back (pair<GlobalOrdinal, GlobalOrdinal> (row, *j));
        sendSizes[numSends]++;
      }
    }
    if (IdsAndRows.size() > 0) {
      numSends++; // one last increment, to make it a count instead of an index
    }
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      as<typename Array<int>::size_type> (numSends) != sendIDs.size (),
      std::logic_error, ": internal logic error. Contact Tpetra team.");

    // don't need this data anymore
    nonlocals_.clear();

    //////////////////////////////////////////////////////////////////////////////////////
    // TRANSMIT SIZE INFO BETWEEN SENDERS AND RECEIVERS
    //////////////////////////////////////////////////////////////////////////////////////

    // Array of pending nonblocking communication requests.  It's OK
    // to mix nonblocking send and receive requests in the same
    // waitAll() call.
    Array<RCP<Teuchos::CommRequest<int> > > requests;

    // perform non-blocking sends: send sizes to our recipients
    for (size_t s = 0; s < numSends ; ++s) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (isend<int, size_t> (*comm,
                                              rcp (&sendSizes[s], false),
                                              sendIDs[s]));
    }
    // perform non-blocking receives: receive sizes from our senders
    Array<size_t> recvSizes (numRecvs);
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      requests.push_back (ireceive (*comm, rcp (&recvSizes[r], false), recvIDs[r]));
    }
    // Wait on all the nonblocking sends and receives.
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    // This doesn't necessarily deallocate the array.
    requests.resize (0);

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW SEND/RECEIVE ALL ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // from the size info, build the ArrayViews into IJSendBuffer
    Array<ArrayView<pair<GO,GO> > > sendBuffers(numSends,null);
    {
      size_t cur = 0;
      for (size_t s=0; s<numSends; ++s) {
        sendBuffers[s] = IJSendBuffer(cur,sendSizes[s]);
        cur += sendSizes[s];
      }
    }
    // perform non-blocking sends
    for (size_t s=0; s < numSends ; ++s)
    {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpSendBuf =
        arcp (sendBuffers[s].getRawPtr(), 0, sendBuffers[s].size(), false);
      requests.push_back (isend<int, pair<GO,GO> > (*comm, tmpSendBuf, sendIDs[s]));
    }
    // calculate amount of storage needed for receives
    // setup pointers for the receives as well
    size_t totalRecvSize = std::accumulate (recvSizes.begin(), recvSizes.end(), 0);
    Array<pair<GO,GO> > IJRecvBuffer (totalRecvSize);
    // from the size info, build the ArrayViews into IJRecvBuffer
    Array<ArrayView<pair<GO,GO> > > recvBuffers (numRecvs, null);
    {
      size_t cur = 0;
      for (size_t r=0; r<numRecvs; ++r) {
        recvBuffers[r] = IJRecvBuffer(cur,recvSizes[r]);
        cur += recvSizes[r];
      }
    }
    // perform non-blocking recvs
    for (size_t r = 0; r < numRecvs; ++r) {
      // We're using a nonowning RCP because all communication
      // will be local to this method and the scope of our data
      ArrayRCP<pair<GO,GO> > tmpRecvBuf =
        arcp (recvBuffers[r].getRawPtr(), 0, recvBuffers[r].size(), false);
      requests.push_back (ireceive (*comm, tmpRecvBuf, recvIDs[r]));
    }
    // perform waits
    if (! requests.empty()) {
      waitAll (*comm, requests());
    }
#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier (*comm);
#endif // HAVE_TPETRA_DEBUG

    ////////////////////////////////////////////////////////////////////////////////////
    // NOW PROCESS THE RECEIVED ROW DATA
    ////////////////////////////////////////////////////////////////////////////////////
    // TODO: instead of adding one entry at a time, add one row at a time.
    //       this requires resorting; they arrived sorted by sending node,
    //       so that entries could be non-contiguous if we received
    //       multiple entries for a particular row from different processors.
    //       it also requires restoring the data, which may make it not worth the trouble.
    for (typename Array<pair<GO,GO> >::const_iterator ij = IJRecvBuffer.begin();
         ij != IJRecvBuffer.end(); ++ij)
    {
      insertGlobalIndicesFiltered (ij->first, tuple<GO> (ij->second));
    }
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  resumeFill (const RCP<ParameterList> &params)
  {
    const char tfecfFuncName[] = "resumeFill";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(! hasRowInfo(), std::runtime_error,
      ": Sorry, you cannot resume fill of the CrsGraph, since the graph's row "
      "information was deleted in fillComplete().");

#ifdef HAVE_TPETRA_DEBUG
    Teuchos::barrier( *rowMap_->getComm() );
#endif // HAVE_TPETRA_DEBUG
    clearGlobalConstants();
    lclGraph_ = null;
    if (params != null) this->setParameterList (params);
    lowerTriangular_  = false;
    upperTriangular_  = false;
    // either still sorted/merged or initially sorted/merged
    indicesAreSorted_ = true;
    noRedundancies_ = true;
    fillComplete_ = false;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      ! isFillActive() || isFillComplete(), std::logic_error,
      "::resumeFill(): At end of method, either fill is not active or fill is "
      "complete.  This violates stated post-conditions.  Please report this bug "
      "to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::fillComplete(const RCP<ParameterList> &params)
  {
    fillComplete(rowMap_,rowMap_,params);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  fillComplete (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &domainMap,
                const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &rangeMap,
                const RCP<ParameterList> &params)
  {
#ifdef HAVE_TPETRA_DEBUG
    rowMap_->getComm ()->barrier ();
#endif // HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "fillComplete()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( ! isFillActive() || isFillComplete(),
      std::runtime_error, ": Graph fill state must be active (isFillActive() "
      "must be true) before calling fillComplete().");

    const int numProcs = getComm ()->getSize ();

    // allocate if unallocated
    if (! indicesAreAllocated()) {
      if (hasColMap ()) {
        // We have a column Map, so use local indices.
        allocateIndices (LocalIndices);
      } else {
        // We don't have a column Map, so use global indices.
        allocateIndices (GlobalIndices);
      }
    }

    // If true, the caller promises that no process did nonlocal
    // changes since the last call to fillComplete.
    bool assertNoNonlocalInserts = false;
    if (! params.is_null ()) {
      assertNoNonlocalInserts =
        params->get<bool> ("No Nonlocal Changes", assertNoNonlocalInserts);
    }
    // We also don't need to do global assembly if there is only one
    // process in the communicator.
    const bool mayNeedGlobalAssemble = ! assertNoNonlocalInserts && numProcs > 1;
    if (mayNeedGlobalAssemble) {
      // This first checks if we need to do global assembly.
      // The check costs a single all-reduce.
      globalAssemble();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
        numProcs > 1 && nonlocals_.size() > 0, std::runtime_error,
        ":" << std::endl << "The graph's communicator contains only one "
        "process, but there are nonlocal entries.  " << std::endl <<
        "This probably means that invalid entries were added to the graph.");
    }
    // set domain/range map: may clear the import/export objects
    setDomainRangeMaps(domainMap,rangeMap);
    // make column map
    if (! hasColMap()) {
      makeColMap();
    }

    // Make indices local, if they aren't already.
    // The method doesn't do any work if the indices are already local.
    makeIndicesLocal ();

    if (! isSorted()) {
      // If this process has no indices, then CrsGraph considers it
      // already trivially sorted.  Thus, this method need not be
      // called on all processes in the row Map's communicator.
      sortAllIndices();
    }

    if (! isMerged()) {
      mergeAllIndices();
    }

    makeImportExport(); // Make Import and Export objects
    computeGlobalConstants();
    // fill local objects
    fillLocalGraph(params);
    //
    fillComplete_ = true;
#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == true || isFillComplete() == false, std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    //
    checkInternalState();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::expertStaticFillComplete(const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > & domainMap,
                                                                                       const RCP<const Map<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > & rangeMap,
                                                                                       const RCP<const Import<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &importer,
                                                                                       const RCP<const Export<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> > > &exporter,
                                                                                       const RCP<ParameterList> &params)
  {
#ifdef HAVE_TPETRA_DEBUG
    rowMap_->getComm ()->barrier ();
#endif // HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "expertStaticFillComplete()";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillComplete() == true || hasColMap() == false,
                                           std::runtime_error, ": fillComplete cannot have already been called and a ColMap is required.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( getNodeNumRows() > 0 && rowPtrs_==Teuchos::null,
                                           std::runtime_error, ": a matrix will getNodeNumRows()>0 requires rowptr to be set.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( domainMap == Teuchos::null || rangeMap == Teuchos::null,
                                           std::runtime_error, ": requires a non-null domainMap & rangeMap.");

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( pftype_ !=StaticProfile,
                                           std::runtime_error, ": requires StaticProfile.");

    // Note: We don't need to do the following things which are normally done in fillComplete:
    // allocateIndices, globalAssemble, makeColMap, makeIndicesLocal, sortAllIndices, mergeAllIndices

    // Note: Need to do this so computeGlobalConstants & fillLocalGraph work
    nodeNumEntries_ = nodeNumAllocated_ = rowPtrs_[getNodeNumRows()];

    // Constants from allocateIndices
    numAllocForAllRows_  = 0;
    numAllocPerRow_      = null;
    indicesAreAllocated_ = true;

    // Constants from makeIndicesLocal
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;

    // set domain/range map: may clear the import/export objects
    setDomainRangeMaps(domainMap,rangeMap);

    // Presume the user sorted and merged the arrays first
    indicesAreSorted_ = true;
    noRedundancies_ = true;

    // makeImportExport won't create a new importer/exporter if I set one here first.
    importer_=Teuchos::null;
    exporter_=Teuchos::null;
    if(importer != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!importer->getSourceMap()->isSameAs(*getDomainMap()) || !importer->getTargetMap()->isSameAs(*getColMap()),
                                            std::invalid_argument,": importer does not match matrix maps.");
      importer_ = importer;

    }
    if(exporter != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(!exporter->getSourceMap()->isSameAs(*getRowMap()) || !exporter->getTargetMap()->isSameAs(*getRangeMap()),
                                            std::invalid_argument,": exporter does not match matrix maps.");
      exporter_ = exporter;
    }
    makeImportExport();

    // Compute the constants
    computeGlobalConstants();

    // Since we have a StaticProfile, fillLocalGraph will do the right thing...
    fillLocalGraph(params);
    fillComplete_ = true;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC( isFillActive() == true || isFillComplete() == false, std::logic_error, ": Violated stated post-conditions. Please contact Tpetra team.");
#endif
    checkInternalState();
  }



  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,
                GlobalOrdinal,
                Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
                typename KokkosClassic::DefaultKernels<void, LocalOrdinal, Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  fillLocalGraph (const RCP<ParameterList> &params)
  {
    const size_t numRows = getNodeNumRows();
    ArrayRCP<size_t> ptrs;
    ArrayRCP<LocalOrdinal> inds;
    t_LocalOrdinal_1D k_inds;
    t_RowPtrs k_ptrs;

    bool requestOptimizedStorage = true;
    if (params != null && params->get("Optimize Storage",true) == false) requestOptimizedStorage = false;
    if (getProfileType() == DynamicProfile) {
      // 2d -> 1d packed
      k_ptrs = t_RowPtrs("Tpetra::CrsGraph::RowPtrs",numRowEntries_.size()+1);
      ptrs = Teuchos::arcp(k_ptrs.ptr_on_device(), 0, k_ptrs.dimension_0(),
                           Kokkos::Compat::deallocator(k_ptrs), false);

      // hack until we get parallel_scan in kokkos
      typename t_RowPtrs::HostMirror h_k_ptrs = Kokkos::create_mirror_view(k_ptrs);
      Kokkos::deep_copy(h_k_ptrs,k_ptrs);
      for(int i = 0; i < numRowEntries_.size(); i++) {
        h_k_ptrs(i+1) = h_k_ptrs(i)+numRowEntries_[i];
      }
      Kokkos::deep_copy(k_ptrs,h_k_ptrs);

      k_inds = t_LocalOrdinal_1D("Tpetra::CrsGraph::lclInds1D_",h_k_ptrs[h_k_ptrs.dimension_0()-1]);
      inds = Teuchos::arcp(k_inds.ptr_on_device(), 0, k_inds.dimension_0(),
                           Kokkos::Compat::deallocator(k_inds), false);
      typename t_LocalOrdinal_1D::HostMirror h_inds = Kokkos::create_mirror_view(k_inds);
      for (size_t row=0; row < numRows; ++row) {
        const size_t numentrs = numRowEntries_[row];
        std::copy (lclInds2D_[row].begin(),
                   lclInds2D_[row].begin() + numentrs,
                   h_inds.ptr_on_device() + h_k_ptrs[row]);
      }
      Kokkos::deep_copy(k_inds,h_inds);
    }

    else if (getProfileType() == StaticProfile) {
      // 1d non-packed -> 1d packed
      if (nodeNumEntries_ != nodeNumAllocated_) {
        k_ptrs = t_RowPtrs("Tpetra::CrsGraph::RowPtrs",numRowEntries_.size()+1);
        ptrs = Teuchos::arcp(k_ptrs.ptr_on_device(), 0, k_ptrs.dimension_0(),
                             Kokkos::Compat::deallocator(k_ptrs), false);
        {
          // hack until we get parallel_scan in kokkos
          typename t_RowPtrs::HostMirror h_k_ptrs = Kokkos::create_mirror_view(k_ptrs);
          Kokkos::deep_copy(h_k_ptrs,k_ptrs);
          for(int i = 0; i < numRowEntries_.size(); i++) {
            h_k_ptrs(i+1) = h_k_ptrs(i)+numRowEntries_[i];
          }
          Kokkos::deep_copy(k_ptrs,h_k_ptrs);
        }

        k_inds = t_LocalOrdinal_1D("Tpetra::CrsGraph::lclInds1D_",*(ptrs.end()-1));
        inds = Teuchos::arcp(k_inds.ptr_on_device(), 0, k_inds.dimension_0(),
                             Kokkos::Compat::deallocator(k_inds), false);
        if (k_rowPtrs_.dimension_0()==0) {
          t_RowPtrs tmp_unpacked_ptrs = t_RowPtrs("Tpetra::CrsGraph::RowPtrs",numRowEntries_.size()+1);
          for(int i = 0; i < rowPtrs_.size(); i++) {
            tmp_unpacked_ptrs(i) = rowPtrs_[i];
          }
          k_rowPtrs_ = tmp_unpacked_ptrs;
        }
        t_LocalOrdinal_1D tmp_unpacked_inds;
        if (k_lclInds1D_.dimension_0()==0) {
          tmp_unpacked_inds = t_LocalOrdinal_1D("Tpetra::CrsGraph::t_LocalOridnal",lclInds1D_.size());
          for(int i = 0; i < lclInds1D_.size(); i++) {
            tmp_unpacked_inds(i) = lclInds1D_[i];
          }
        } else
          tmp_unpacked_inds = k_lclInds1D_;
        {
          pack_functor<t_LocalOrdinal_1D, t_RowPtrs>
            f(k_inds,tmp_unpacked_inds,k_ptrs,k_rowPtrs_);
          Kokkos::parallel_for(numRows,f);
        }

      }
      else {
        inds = lclInds1D_;
        ptrs = rowPtrs_;
        k_ptrs = k_rowPtrs_;
        k_inds = k_lclInds1D_;
      }
    }
    // can we ditch the old allocations for the packed one?
    if ( requestOptimizedStorage ) {
      lclInds2D_ = null;
      numRowEntries_ = null;
      // keep the new stuff
      lclInds1D_ = inds;
      rowPtrs_ = ptrs;
      nodeNumAllocated_ = nodeNumEntries_;
      pftype_ = StaticProfile;
      k_lclInds1D_ = k_inds;
      k_rowPtrs_   = k_ptrs;
    }
    // build the local graph, hand over the indices
    RCP<ParameterList> lclparams;
    if (params == null) lclparams = parameterList();
    else                lclparams = sublist(params,"Local Graph");
    lclGraph_ = rcp( new local_graph_type( getRowMap()->getNodeNumElements(), getColMap()->getNodeNumElements(), getRowMap()->getNode(), lclparams ) );
    lclGraph_->setStructure(ptrs,inds);
    k_lclGraph_ = LocalStaticCrsGraphType(k_inds,k_ptrs);
    ptrs = null;
    inds = null;
    // finalize local graph
    //
    // TODO (mfh 13 Mar 2014) getNodeNumDiags(), isUpperTriangular(),
    // and isLowerTriangular() depend on computeGlobalConstants(), in
    // particular the part where it looks at the local matrix.  You
    // have to use global indices to determine which entries are
    // diagonal, or above or below the diagonal.  However, lower or
    // upper triangularness is a local property.
    Teuchos::EDiag diag = ( getNodeNumDiags() < getNodeNumRows() ? Teuchos::UNIT_DIAG : Teuchos::NON_UNIT_DIAG );
    Teuchos::EUplo uplo = Teuchos::UNDEF_TRI;
    if      (isUpperTriangular()) uplo = Teuchos::UPPER_TRI;
    else if (isLowerTriangular()) uplo = Teuchos::LOWER_TRI;
    LocalMatOps::finalizeGraph(uplo,diag,*lclGraph_,params);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
    typename KokkosClassic::DefaultKernels<
      void,
      LocalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  replaceColMap (const Teuchos::RCP<const map_type>& newColMap)
  {
    // NOTE: This safety check matches the code, but not the documentation of Crsgraph
    const char tfecfFuncName[] = "replaceColMap";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(isLocallyIndexed() || isGloballyIndexed(),  std::runtime_error, " requires matching maps and non-static graph.");
    colMap_ = newColMap;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
    typename KokkosClassic::DefaultKernels<
      void,
      LocalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  replaceDomainMapAndImporter (const Teuchos::RCP<const map_type>& newDomainMap,
                               const Teuchos::RCP<const import_type>& newImporter)
  {
    const char prefix[] = "Tpetra::CrsGraph::replaceDomainMapAndImporter: ";
    TEUCHOS_TEST_FOR_EXCEPTION(
      colMap_.is_null (), std::invalid_argument, prefix << "You may not call "
      "this method unless the graph already has a column Map.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      newDomainMap.is_null (), std::invalid_argument,
      prefix << "The new domain Map must be nonnull.");

#ifdef HAVE_TPETRA_DEBUG
    if (newImporter.is_null ()) {
      // It's not a good idea to put expensive operations in a macro
      // clause, even if they are side effect - free, because macros
      // don't promise that they won't evaluate their arguments more
      // than once.  It's polite for them to do so, but not required.
      const bool colSameAsDom = colMap_->isSameAs (*newDomainMap);
      TEUCHOS_TEST_FOR_EXCEPTION(
        colSameAsDom, std::invalid_argument, "If the new Import is null, "
        "then the new domain Map must be the same as the current column Map.");
    }
    else {
      const bool colSameAsTgt =
        colMap_->isSameAs (* (newImporter->getTargetMap ()));
      const bool newDomSameAsSrc =
        newDomainMap->isSameAs (* (newImporter->getSourceMap ()));
      TEUCHOS_TEST_FOR_EXCEPTION(
        ! colSameAsTgt || ! newDomSameAsSrc, std::invalid_argument, "If the "
        "new Import is nonnull, then the current column Map must be the same "
        "as the new Import's target Map, and the new domain Map must be the "
        "same as the new Import's source Map.");
    }
#endif // HAVE_TPETRA_DEBUG

    domainMap_ = newDomainMap;
    importer_ = Teuchos::rcp_const_cast<import_type> (newImporter);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  const RCP<const typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps::template graph<LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::graph_type>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getLocalGraph() const
  {
    return lclGraph_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  Kokkos::StaticCrsGraph<LocalOrdinal, Kokkos::LayoutLeft, typename Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>::device_type,size_t>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getLocalGraph_Kokkos() const
  {
    return k_lclGraph_;
  }

  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  const RCP<typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps::template graph<LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::graph_type>
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::getLocalGraphNonConst()
  {
    return lclGraph_;
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<
    LocalOrdinal,
    GlobalOrdinal,
    Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
    typename KokkosClassic::DefaultKernels<
      void,
      LocalOrdinal,
      Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  computeGlobalConstants ()
  {
    using Teuchos::as;
    using Teuchos::outArg;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;

#ifdef HAVE_TPETRA_DEBUG
    TEUCHOS_TEST_FOR_EXCEPTION(! hasColMap(), std::logic_error, "Tpetra::"
      "CrsGraph::computeGlobalConstants: At this point, the graph should have "
      "a column Map, but it does not.  Please report this bug to the Tpetra "
      "developers.");
#endif // HAVE_TPETRA_DEBUG

    // If necessary, (re)compute the local constants: nodeNumDiags_,
    // lowerTriangular_, upperTriangular_, and nodeMaxNumRowEntries_.
    if (! haveLocalConstants_) {
      // We have actually already computed nodeNumEntries_.
      // nodeNumEntries_ gets updated by methods that insert or remove
      // indices (including setAllIndices and
      // expertStaticFillComplete).  Before fillComplete, its count
      // may include duplicate column indices in the same row.
      // However, mergeRowIndices and mergeRowIndicesAndValues both
      // subtract off merged indices in each row from the total count.
      // Thus, nodeNumEntries_ _should_ be accurate at this point,
      // meaning that we don't have to re-count it here.

      // Reset local properties
      upperTriangular_ = true;
      lowerTriangular_ = true;
      nodeMaxNumRowEntries_ = 0;
      nodeNumDiags_         = 0;

      // At this point, we know that we have both a row Map and a column Map.
      const Map<LO,GO,Node>& rowMap = *rowMap_;
      const Map<LO,GO,Node>& colMap = *colMap_;

      // Go through all the entries of the graph.  Count the number of
      // diagonal elements we encounter, and figure out whether the
      // graph is lower or upper triangular.  Diagonal elements are
      // determined using global indices, with respect to the whole
      // graph.  However, lower or upper triangularity is a local
      // property, and is determined using local indices.
      //
      // At this point, indices have already been sorted in each row.
      // That makes finding out whether the graph is lower / upper
      // triangular easier.
      if (indicesAreAllocated () && nodeNumAllocated_ > 0) {
        const size_t numLocalRows = getNodeNumRows ();
        for (size_t localRow = 0; localRow < numLocalRows; ++localRow) {
          const GO globalRow = rowMap.getGlobalElement (localRow);
          // Find the local (column) index for the diagonal entry.
          // This process might not necessarily own _any_ entries in
          // the current row.  If it doesn't, skip this row.  It won't
          // affect any of the attributes (nodeNumDiagons_,
          // upperTriangular_, lowerTriangular_, or
          // nodeMaxNumRowEntries_) which this loop sets.
          const LO rlcid = colMap.getLocalElement (globalRow);
          if (rlcid != Teuchos::OrdinalTraits<LO>::invalid ()) {
            // This process owns one or more entries in the current row.
            RowInfo rowInfo = getRowInfo (localRow);
            ArrayView<const LO> rview = getLocalView (rowInfo);
            typename ArrayView<const LO>::iterator beg, end, cur;
            beg = rview.begin();
            end = beg + rowInfo.numEntries;
            if (beg != end) {
              for (cur = beg; cur != end; ++cur) {
                // is this the diagonal?
                if (rlcid == *cur) ++nodeNumDiags_;
              }
              // Local column indices are sorted in each row.  That means
              // the smallest column index in this row (on this process)
              // is *beg, and the largest column index in this row (on
              // this process) is *(end - 1).  We know that end - 1 is
              // valid because beg != end.
              const size_t smallestCol = as<size_t> (*beg);
              const size_t largestCol = as<size_t> (*(end - 1));

              if (smallestCol < localRow) {
                upperTriangular_ = false;
              }
              if (localRow < largestCol) {
                lowerTriangular_ = false;
              }
            }
            // Update the max number of entries over all rows.
            nodeMaxNumRowEntries_ = std::max (nodeMaxNumRowEntries_, rowInfo.numEntries);
          }
        }
      }
      haveLocalConstants_ = true;
    } // if my process doesn't have local constants

    // Compute global constants from local constants.  Processes that
    // already have local constants still participate in the
    // all-reduces, using their previously computed values.
    if (haveGlobalConstants_ == false) {
      // Promote all the nodeNum* and nodeMaxNum* quantities from
      // size_t to global_size_t, when doing the all-reduces for
      // globalNum* / globalMaxNum* results.
      //
      // FIXME (mfh 07 May 2013) Unfortunately, we either have to do
      // this in two all-reduces (one for the sum and the other for
      // the max), or use a custom MPI_Op that combines the sum and
      // the max.  The latter might even be slower than two
      // all-reduces on modern network hardware.  It would also be a
      // good idea to use nonblocking all-reduces (MPI 3), so that we
      // don't have to wait around for the first one to finish before
      // starting the second one.
      GST lcl[2], gbl[2];
      lcl[0] = as<GST> (nodeNumEntries_);
      lcl[1] = as<GST> (nodeNumDiags_);
      reduceAll<int,GST> (*getComm (), Teuchos::REDUCE_SUM,
                          2, lcl, gbl);
      globalNumEntries_ = gbl[0];
      globalNumDiags_   = gbl[1];
      reduceAll<int,GST> (*getComm (), Teuchos::REDUCE_MAX,
                          as<GST> (nodeMaxNumRowEntries_),
                          outArg (globalMaxNumRowEntries_));
      haveGlobalConstants_ = true;
    }
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal, GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<void, LocalOrdinal,
                                                  Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  makeIndicesLocal ()
  {
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
#ifdef HAVE_TPETRA_DEBUG
    const char tfecfFuncName[] = "makeIndicesLocal";
#endif // HAVE_TPETRA_DEBUG

    // Transform indices to local index space
    const size_t nlrs = getNodeNumRows ();
    if (isGloballyIndexed () && nlrs > 0) {
      // allocate data for local indices
      if (getProfileType() == StaticProfile) {
        // If GO and LO are the same size, we can reuse the existing
        // array of 1-D index storage to convert column indices from
        // GO to LO.  Otherwise, we'll just allocate a new buffer.
        if (nodeNumAllocated_ && Kokkos::Impl::is_same<LO,GO>::value) {
          lclInds1D_ = arcp_reinterpret_cast<LO> (gblInds1D_).persistingView (0, nodeNumAllocated_);
          k_lclInds1D_ = Kokkos::Impl::if_c<Kokkos::Impl::is_same<LO,GO>::value,
            t_GlobalOrdinal_1D,
            t_LocalOrdinal_1D >::select (k_gblInds1D_, k_lclInds1D_);
        }
        else {
          k_lclInds1D_ = t_LocalOrdinal_1D ("Tpetra::CrsGraph::lclInds1D_", * (rowPtrs_.end () - 1));
          lclInds1D_ = arcp (k_lclInds1D_.ptr_on_device (), 0,
                             k_lclInds1D_.dimension_0 (),
                             Kokkos::Compat::deallocator (k_lclInds1D_),
                             false);
        }

        for (size_t r = 0; r < nlrs; ++r) {
          const size_t offset   = rowPtrs_[r];
          const size_t numentry = numRowEntries_[r];
          for (size_t j = 0; j < numentry; ++j) {
            const GO gid = gblInds1D_[offset + j];
            const LO lid = colMap_->getLocalElement (gid);
            lclInds1D_[offset + j] = lid;
#ifdef HAVE_TPETRA_DEBUG
            TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
              lclInds1D_[offset + j] == Teuchos::OrdinalTraits<LO>::invalid(),
              std::logic_error,
              ": In local row r=" << r << ", global column " << gid << " is "
              "not in the column Map.  This should never happen.  Please "
              "report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
          }
        }
        // We've converted column indices from global to local, so we
        // can deallocate the global column indices (which we know are
        // in 1-D storage, because the graph has static profile).
        gblInds1D_ = null;
      }
      else {  // the graph has dynamic profile (2-D index storage)
        lclInds2D_ = arcp<Array<LO> > (nlrs);
        for (size_t r = 0; r < getNodeNumRows (); ++r) {
          if (! gblInds2D_[r].empty ()) {
            const GO* const ginds = gblInds2D_[r].getRawPtr ();
            const size_t rna = gblInds2D_[r].size ();
            const size_t numentry = numRowEntries_[r];
            lclInds2D_[r].resize (rna);
            LO* const linds = lclInds2D_[r].getRawPtr ();
            for (size_t j = 0; j < numentry; ++j) {
              GO gid = ginds[j];
              LO lid = colMap_->getLocalElement (gid);
              linds[j] = lid;
#ifdef HAVE_TPETRA_DEBUG
              TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
                linds[j] == Teuchos::OrdinalTraits<LO>::invalid(),
                std::logic_error,
                ": Global column ginds[j=" << j << "]=" << ginds[j]
                << " of local row r=" << r << " is not in the column Map.  "
                "This should never happen.  Please report this bug to the "
                "Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
            }
          }
        }
        gblInds2D_ = null;
      }
    } // globallyIndexed() && nlrs > 0

    k_lclGraph_ = LocalStaticCrsGraphType (k_lclInds1D_, k_rowPtrs_);
    indicesAreLocal_  = true;
    indicesAreGlobal_ = false;
    checkInternalState ();
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,
           GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<
             void,
             LocalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  computeIndexState ()
  {
    // FIXME (mfh 03 Mar 2013) It's not clear to me that we need to do
    // an all-reduce.
    //
    // This method is _only_ called in makeIndicesLocal() and
    // makeColMap().  makeIndicesLocal() calls makeColMap(), which is
    // collective, so both methods must be called collectively.
    //
    // There are only two methods that modify indicesAreLocal_:
    // allocateIndices(), and makeIndicesLocal().  makeIndicesLocal()
    // always has the eponymous effect.  It must be called
    // collectively, since it calls makeColMap(), which must be called
    // collectively.
    //
    // allocateIndices(), on the other hand, could perhaps be called
    // with lg=LocalIndices on one process, but with lg=GlobalIndices
    // on another.  However, I will argue that CrsGraph cannot reach
    // this state if used correctly.
    //
    // allocateIndices() specifically forbids an lg argument which
    // does not correspond to the state of the graph.  The graph
    // starts out locally indexed if its constructor was given a
    // column Map; otherwise, it starts out neither locally nor
    // globally indexed, and only gains one of these identities on the
    // calling process ("locally") once the user inserts an entry.
    // (This is the classic Epetra way to tell if a graph is empty.)
    // It is an error to call different constructors for the same
    // CrsGraph instance on different processes.  Thus, the initial
    // local-or-global state before any insertions on any processes
    // must be the same.
    //
    // Before inserting any entries, indicesAreLocal_ and
    // indicesAreGlobal_ are both locally false.  They may be modified
    // locally by calls to insertGlobalIndices() or
    // insertLocalIndices().  These two methods only call
    // allocateIndices() if no insertion method has yet been called on
    // the graph.  Furthermore, these methods only allow indices to
    // have the state matching their name.  insertLocalIndices()
    // explicitly requires that the graph has a column Map.
    // insertGlobalIndices() requires that the graph not be locally
    // indexed, which currently means that the graph was not
    // constructed with a column Map and that fillComplete() (which is
    // collective) has not yet been called.
    //
    // Thus, before calling fillComplete() for the first time, it is
    // possible that on some (but not necessarily all) processes,
    // indicesAreLocal_ == false && indicesAreGlobal_ == false.
    // However, on all processes p on which any one of these Booleans
    // are true, the two Booleans must have the same values, for the
    // reasons discussed in the previous paragraph.
    //
    // fillComplete() makes the column Map first (if it doesn't
    // already exist) before making indices local, so there is a point
    // in fillComplete() at which the graph has a column Map and is
    // globally indexed.  However, makeIndicesLocal() fixes this state
    // right away.  This intermediate state would never be exposed to
    // users.  fillComplete() must be called collectively.
    //
    // resumeFill() does _not_ currently change the global vs. local
    // indices state of the graph.  If we were to give users the
    // option to do this in the future (e.g., they want to insert
    // column indices not in the column Map, so that we would have to
    // convert all the column indices back to global first and get rid
    // of the existing column Map), then that would be a collective
    // decision in any case.
    //
    // The only part of makeIndicesLocal() that is not a local
    // operation is the call to makeColMap().  Everything else is
    // local.  Thus, it suffices in that method to check the local
    // state.  Furthermore, makeIndicesLocal() always makes indices
    // local, so there is no need to check at the end of the method.
    // One would only call makeColMap() if the graph does not have a
    // column Map.  In that case, the graph must be globally indexed
    // anyway.
    int myIndices[2] = {0,0};
    if (indicesAreLocal_)  myIndices[0] = 1;
    if (indicesAreGlobal_) myIndices[1] = 1;
    int allIndices[2];
    Teuchos::reduceAll<int, int> (* (getComm()), Teuchos::REDUCE_MAX,
                                  2, myIndices, allIndices);
    // If indices are (local, global) on one process, they should be
    // (local, global) on all processes.
    indicesAreLocal_  = (allIndices[0]==1);
    indicesAreGlobal_ = (allIndices[1]==1);
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,
           GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<
             void,
             LocalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  sortAllIndices ()
  {
    TEUCHOS_TEST_FOR_EXCEPT(isGloballyIndexed()==true);   // this should be called only after makeIndicesLocal()
    if (isSorted () == false) {
      // FIXME (mfh 06 Mar 2014) This would be a good place for a
      // thread-parallel kernel.
      for (size_t row = 0; row < getNodeNumRows (); ++row) {
        RowInfo rowInfo = getRowInfo (row);
        sortRowIndices (rowInfo);
      }
    }
    indicesAreSorted_ = true; // we just sorted every row
  }


  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,
           GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<
             void,
             LocalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  makeColMap ()
  {
    using std::endl;
    using Teuchos::REDUCE_MAX;
    using Teuchos::reduceAll;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "makeColMap";

    if (hasColMap ()) { // The graph already has a column Map.
      // FIXME (mfh 26 Feb 2013): This currently prevents structure
      // changes that affect the column Map.
      return;
    }

    computeIndexState ();
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isLocallyIndexed (), std::runtime_error,
      ": The graph is locally indexed.  Calling makeColMap() to make the "
      "column Map requires that the graph be globally indexed.");

    // After the calling process is done going through all of the rows
    // it owns, myColumns will contain the list of indices owned by
    // this process in the column Map.
    Array<GO> myColumns;

    // If we reach this point, we don't have a column Map yet, so the
    // graph can't be locally indexed.  Thus, isGloballyIndexed() ==
    // false means that the graph is empty on this process, so
    // myColumns will be left empty.
    if (isGloballyIndexed ()) {
      // Go through all the rows, finding the populated column indices.
      //
      // Our final list of indices for the column Map constructor will
      // have the following properties (all of which are with respect
      // to the calling process):
      //
      // 1. Indices in the domain Map go first.
      // 2. Indices not in the domain Map follow, ordered first
      //    contiguously by their owning process rank (in the domain
      //    Map), then in increasing order within that.
      // 3. No duplicate indices.
      //
      // This imitates the ordering used by Aztec(OO) and Epetra.
      // Storing indices owned by the same process (in the domain Map)
      // contiguously permits the use of contiguous send and receive
      // buffers.
      //
      // We begin by partitioning the column indices into "local" GIDs
      // (owned by the domain Map) and "remote" GIDs (not owned by the
      // domain Map).  We use the same order for local GIDs as the
      // domain Map, so we track them in place in their array.  We use
      // an std::set (RemoteGIDSet) to keep track of remote GIDs, so
      // that we don't have to merge duplicates later.
      const LO LINV = Teuchos::OrdinalTraits<LO>::invalid ();
      size_t numLocalColGIDs = 0, numRemoteColGIDs = 0;

      // GIDisLocal[lid] == 0 if and only if local index lid in the
      // domain Map is remote (not local).
      Array<char> GIDisLocal (domainMap_->getNodeNumElements (), 0);
      std::set<GO> RemoteGIDSet;
      const size_t myNumRows = getNodeNumRows ();
      for (size_t r = 0; r < myNumRows; ++r) {
        RowInfo rowinfo = getRowInfo (r);
        if (rowinfo.numEntries > 0) {
          // FIXME (mfh 03 Mar 2013) It's a bit puzzling to me why the
          // ArrayView that getGlobalView() returns doesn't return
          // rowinfo.numEntries entries.
          ArrayView<const GO> rowGids = getGlobalView (rowinfo);
          rowGids = rowGids (0, rowinfo.numEntries);

          for (size_t k = 0; k < rowinfo.numEntries; ++k) {
            const GO gid = rowGids[k];
            const LO lid = domainMap_->getLocalElement (gid);
            if (lid != LINV) {
              const char alreadyFound = GIDisLocal[lid];
              if (alreadyFound == 0) {
                GIDisLocal[lid] = 1;
                ++numLocalColGIDs;
              }
            }
            else {
              const bool notAlreadyFound = RemoteGIDSet.insert (gid).second;
              if (notAlreadyFound) { // gid did not exist in the set before
                ++numRemoteColGIDs;
              }
            }
          } // for each entry k in row r
        } // if row r contains > 0 entries
      } // for each locally owned row r

      // Possible short-circuit for serial scenario:
      //
      // If all domain GIDs are present as column indices, then set
      // ColMap=DomainMap.  By construction, LocalGIDs is a subset of
      // DomainGIDs.
      //
      // If we have
      //   * Number of remote GIDs is 0, so that ColGIDs == LocalGIDs,
      // and
      //   * Number of local GIDs is number of domain GIDs
      // then
      //   * LocalGIDs \subset DomainGIDs && size(LocalGIDs) ==
      //     size(DomainGIDs) => DomainGIDs == LocalGIDs == ColGIDs
      // on the calling process.
      //
      // We will concern ourselves only with the special case of a
      // serial DomainMap, obviating the need for communication.
      //
      // If
      //   * DomainMap has a serial communicator
      // then we can set the column Map as the domain Map
      // return. Benefit: this graph won't need an Import object
      // later.
      //
      // Note, for a serial domain map, there can be no RemoteGIDs,
      // because there are no remote processes.  Likely explanations
      // for this are:
      //  * user submitted erroneous column indices
      //  * user submitted erroneous domain Map
      if (domainMap_->getComm ()->getSize () == 1) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          numRemoteColGIDs != 0, std::runtime_error,
          ": " << numRemoteColGIDs << " column "
          << (numRemoteColGIDs != 1 ? "indices are" : "index is")
          << " not in the domain Map." << endl
          << "Either these indices are invalid or the domain Map is invalid."
          << endl << "Remember that nonsquare matrices, or matrices where the "
          "row and range Maps are different, require calling the version of "
          "fillComplete that takes the domain and range Maps as input.");
        if (numLocalColGIDs == domainMap_->getNodeNumElements()) {
          colMap_ = domainMap_;
          checkInternalState ();
          return;
        }
      }

      // Populate myColumns with a list of all column GIDs.  Put
      // locally owned (in the domain Map) GIDs at the front: they
      // correspond to "same" and "permuted" entries between the
      // column Map and the domain Map.  Put remote GIDs at the back.
      myColumns.resize (numLocalColGIDs + numRemoteColGIDs);
      // get pointers into myColumns for each part
      ArrayView<GO> LocalColGIDs  = myColumns (0, numLocalColGIDs);
      ArrayView<GO> RemoteColGIDs = myColumns (numLocalColGIDs, numRemoteColGIDs);

      // Copy the remote GIDs into myColumns
      std::copy (RemoteGIDSet.begin(), RemoteGIDSet.end(), RemoteColGIDs.begin());

      // Make a list of process ranks corresponding to the remote GIDs.
      Array<int> RemoteImageIDs (numRemoteColGIDs);
      // Look up the remote process' ranks in the domain Map.
      {
        const LookupStatus stat =
          domainMap_->getRemoteIndexList (RemoteColGIDs, RemoteImageIDs ());
#ifdef HAVE_TPETRA_DEBUG
        // If any process returns IDNotPresent, then at least one of
        // the remote indices was not present in the domain Map.  This
        // means that the Import object cannot be constructed, because
        // of incongruity between the column Map and domain Map.
        // This has two likely causes:
        //   - The user has made a mistake in the column indices
        //   - The user has made a mistake with respect to the domain Map
        const int missingID_lcl = (stat == IDNotPresent ? 1 : 0);
        int missingID_gbl = 0;
        reduceAll<int, int> (*getComm (), REDUCE_MAX, missingID_lcl,
                             outArg (missingID_gbl));
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          missingID_gbl == 1, std::runtime_error,
          ": Some column indices are not in the domain Map." << endl
          << "Either these column indices are invalid or the domain Map is "
          "invalid." << endl << "Likely cause: For a nonsquare matrix, you "
          "must give the domain and range Maps as input to fillComplete.");
#else
        (void) stat; // forestall compiler warning for unused variable
#endif // HAVE_TPETRA_DEBUG
      }
      // Sort incoming remote column indices so that all columns
      // coming from a given remote process are contiguous.  This
      // means the Import's Distributor doesn't need to reorder data.
      sort2 (RemoteImageIDs.begin(), RemoteImageIDs.end(), RemoteColGIDs.begin());

      // Copy the local GIDs into myColumns. Two cases:
      // 1. If the number of Local column GIDs is the same as the
      //    number of Local domain GIDs, we can simply read the domain
      //    GIDs into the front part of ColIndices (see logic above
      //    from the serial short circuit case)
      // 2. We step through the GIDs of the DomainMap, checking to see
      //    if each domain GID is a column GID.  We want to do this to
      //    maintain a consistent ordering of GIDs between the columns
      //    and the domain.

      // FIXME (mfh 03 Mar 2013) It's common that the domain Map is
      // contiguous.  It would be more efficient in that case to avoid
      // calling getNodeElementList(), since that permanently
      // constructs and caches the GID list in the contiguous Map.
      ArrayView<const GO> domainElts = domainMap_->getNodeElementList ();
      const size_t numDomainElts = domainMap_->getNodeNumElements ();
      if (numLocalColGIDs == numDomainElts) {
        // If the number of locally owned GIDs are the same as the
        // number of local domain Map elements, then the local domain
        // Map elements are the same as the locally owned GIDs.
        std::copy (domainElts.begin(), domainElts.end(), LocalColGIDs.begin());
      }
      else {
        // Count the number of locally owned GIDs, both to keep track
        // of the current array index, and as a sanity check.
        size_t numLocalCount = 0;
        for (size_t i = 0; i < numDomainElts; ++i) {
          if (GIDisLocal[i]) {
            LocalColGIDs[numLocalCount++] = domainElts[i];
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
          numLocalCount != numLocalColGIDs, std::logic_error,
          ": numLocalCount = " << numLocalCount << " != numLocalColGIDs = "
          << numLocalColGIDs << ".  This should never happen.  Please report "
          "this bug to the Tpetra developers.");
      }

      // FIXME (mfh 03 Apr 2013) Now would be a good time to use the
      // information we collected above to construct the Import.  In
      // particular, building an Import requires:
      //
      // 1. numSameIDs (length of initial contiguous sequence of GIDs
      //    on this process that are the same in both Maps; this
      //    equals the number of domain Map elements on this process)
      //
      // 2. permuteToLIDs and permuteFromLIDs (both empty in this
      //    case, since there's no permutation going on; the column
      //    Map starts with the domain Map's GIDs, and immediately
      //    after them come the remote GIDs)
      //
      // 3. remoteGIDs (exactly those GIDs that we found out above
      //    were not in the domain Map) and remoteLIDs (which we could
      //    have gotten above by using the three-argument version of
      //    getRemoteIndexList() that computes local indices as well
      //    as process ranks, instead of the two-argument version that
      //    was used above)
      //
      // 4. remotePIDs (which we have from the getRemoteIndexList()
      //    call above)
      //
      // 5. Sorting remotePIDs, and applying that permutation to
      //    remoteGIDs and remoteLIDs (by calling sort3 above instead
      //    of sort2)
      //
      // 6. Everything after the sort3 call in Import::setupExport():
      //    a. Create the Distributor via createFromRecvs(), which
      //       computes exportGIDs and exportPIDs
      //    b. Compute exportLIDs from exportGIDs (by asking the
      //       source Map, in this case the domain Map, to convert
      //       global to local)
      //
      // Steps 1-5 come for free, since we must do that work anyway in
      // order to compute the column Map.  In particular, Step 3 is
      // even more expensive than Step 6a, since it involves both
      // creating and using a new Distributor object.

    } // if the graph is globally indexed

    const global_size_t gstInv =
      Teuchos::OrdinalTraits<global_size_t>::invalid ();
    // FIXME (mfh 05 Mar 2014) Doesn't the index base of a Map have to
    // be the same as the Map's min GID? If the first column is empty
    // (contains no entries), then the column Map's min GID won't
    // necessarily be the same as the domain Map's index base.
    const GO indexBase = domainMap_->getIndexBase ();
    colMap_ = rcp (new map_type (gstInv, myColumns, indexBase,
                                 domainMap_->getComm (),
                                 domainMap_->getNode ()));

    checkInternalState ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::mergeAllIndices()
  {
    TEUCHOS_TEST_FOR_EXCEPT( isGloballyIndexed() ); // call only after makeIndicesLocal()
    TEUCHOS_TEST_FOR_EXCEPT( ! isSorted() ); // call only after sortIndices()
    if (! isMerged ()) {
      for (size_t row=0; row < getNodeNumRows(); ++row) {
        RowInfo rowInfo = getRowInfo(row);
        mergeRowIndices(rowInfo);
      }
      // we just merged every row
      noRedundancies_ = true;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,
           GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<
             void,
             LocalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  makeImportExport ()
  {
    TEUCHOS_TEST_FOR_EXCEPTION(! hasColMap (), std::logic_error, "Tpetra::"
      "CrsGraph::makeImportExport: This method may not be called unless the "
      "graph has a column Map.");
    RCP<ParameterList> params = this->getNonconstParameterList (); // could be null

    // Don't do any checks to see if we need to create the Import, if
    // it exists already.
    //
    // FIXME (mfh 25 Mar 2013) This will become incorrect if we
    // change CrsGraph in the future to allow changing the column
    // Map after fillComplete.  For now, the column Map is fixed
    // after the first fillComplete call.
    if (importer_.is_null ()) {
      // Create the Import instance if necessary.
      if (domainMap_ != colMap_ && (! domainMap_->isSameAs (*colMap_))) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          importer_ = rcp (new import_type (domainMap_, colMap_));
        } else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          importer_ = rcp (new import_type (domainMap_, colMap_, importSublist));
        }
      }
    }

    // Don't do any checks to see if we need to create the Export, if
    // it exists already.
    if (exporter_.is_null ()) {
      // Create the Export instance if necessary.
      if (rangeMap_ != rowMap_ && ! rangeMap_->isSameAs (*rowMap_)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter_ = rcp (new export_type (rowMap_, rangeMap_));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter_ = rcp (new export_type (rowMap_, rangeMap_, exportSublist));
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  std::string
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::description() const
  {
    std::ostringstream oss;
    oss << DistObject<GlobalOrdinal,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::description();
    if (isFillComplete()) {
      oss << "{status = fill complete"
          << ", global rows = " << getGlobalNumRows()
          << ", global cols = " << getGlobalNumCols()
          << ", global num entries = " << getGlobalNumEntries()
          << "}";
    }
    else {
      oss << "{status = fill not complete"
          << ", global rows = " << getGlobalNumRows()
          << "}";
    }
    return oss.str();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    const char tfecfFuncName[] = "describe()";
    using std::endl;
    using std::setw;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    Teuchos::EVerbosityLevel vl = verbLevel;
    if (vl == VERB_DEFAULT) vl = VERB_LOW;
    RCP<const Comm<int> > comm = this->getComm();
    const int myImageID = comm->getRank(),
              numImages = comm->getSize();
    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumRows(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t>(width,Teuchos::as<size_t>(11)) + 2;
    Teuchos::OSTab tab(out);
    //    none: print nothing
    //     low: print O(1) info from node 0
    //  medium: print O(P) info, num entries per node
    //    high: print O(N) info, num entries per row
    // extreme: print O(NNZ) info: print graph indices
    //
    // for medium and higher, print constituent objects at specified verbLevel
    if (vl != VERB_NONE) {
      if (myImageID == 0) out << this->description() << std::endl;
      // O(1) globals, minus what was already printed by description()
      if (isFillComplete() && myImageID == 0) {
        out << "Global number of diagonals = " << globalNumDiags_ << std::endl;
        out << "Global max number of row entries = " << globalMaxNumRowEntries_ << std::endl;
      }
      // constituent objects
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        if (myImageID == 0) out << "\nRow map: " << std::endl;
        rowMap_->describe(out,vl);
        if (colMap_ != null) {
          if (myImageID == 0) out << "\nColumn map: " << std::endl;
          colMap_->describe(out,vl);
        }
        if (domainMap_ != null) {
          if (myImageID == 0) out << "\nDomain map: " << std::endl;
          domainMap_->describe(out,vl);
        }
        if (rangeMap_ != null) {
          if (myImageID == 0) out << "\nRange map: " << std::endl;
          rangeMap_->describe(out,vl);
        }
      }
      // O(P) data
      if (vl == VERB_MEDIUM || vl == VERB_HIGH || vl == VERB_EXTREME) {
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << "Node ID = " << imageCtr << std::endl
                << "Node number of entries = " << nodeNumEntries_ << std::endl
                << "Node number of diagonals = " << nodeNumDiags_ << std::endl
                << "Node max number of entries = " << nodeMaxNumRowEntries_ << std::endl;
            if (indicesAreAllocated()) {
              out << "Node number of allocated entries = " << nodeNumAllocated_ << std::endl;
            }
            else {
              out << "Indices are not allocated." << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
      // O(N) and O(NNZ) data
      if (vl == VERB_HIGH || vl == VERB_EXTREME) {
        TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(hasRowInfo() == false, std::runtime_error, ": reduce verbosity level; graph row information was deleted at fillComplete().");
        for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
          if (myImageID == imageCtr) {
            out << std::setw(width) << "Node ID"
                << std::setw(width) << "Global Row"
                << std::setw(width) << "Num Entries";
            if (vl == VERB_EXTREME) {
              out << "Entries";
            }
            out << std::endl;
            for (size_t r=0; r < getNodeNumRows(); ++r) {
              RowInfo rowinfo = getRowInfo(r);
              GlobalOrdinal gid = rowMap_->getGlobalElement(r);
              out << std::setw(width) << myImageID
                  << std::setw(width) << gid
                  << std::setw(width) << rowinfo.numEntries;
              if (vl == VERB_EXTREME) {
                if (isGloballyIndexed()) {
                  ArrayView<const GlobalOrdinal> rowview = getGlobalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << rowview[j] << " ";
                }
                else if (isLocallyIndexed()) {
                  ArrayView<const LocalOrdinal> rowview = getLocalView(rowinfo);
                  for (size_t j=0; j < rowinfo.numEntries; ++j) out << colMap_->getGlobalElement(rowview[j]) << " ";
                }
              }
              out << std::endl;
            }
          }
          comm->barrier();
          comm->barrier();
          comm->barrier();
        }
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  checkSizes (const SrcDistObject& source)
  {
    (void) source; // forestall "unused variable" compiler warnings

    // It's not clear what kind of compatibility checks on sizes can
    // be performed here.  Epetra_CrsGraph doesn't check any sizes for
    // compatibility.
    return true;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  copyAndPermute (const SrcDistObject& source,
                  size_t numSameIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                  const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs)
  {
    using Teuchos::Array;
    using Teuchos::ArrayView;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "copyAndPermute";
    typedef CrsGraph<LO, GO, Node, LocalMatOps> this_type;
    typedef RowGraph<LO, GO, Node> row_graph_type;

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      permuteToLIDs.size() != permuteFromLIDs.size(), std::runtime_error,
      ": permuteToLIDs and permuteFromLIDs must have the same size.");
    // Make sure that the source object has the right type.  We only
    // actually need it to be a RowGraph, with matching first three
    // template parameters.  If it's a CrsGraph, we can use view mode
    // instead of copy mode to get each row's data.
    //
    // FIXME (mfh 07 Jul 2013) It should not be necessary for any of
    // the template parameters but GO to match.  GO has to match
    // because the graph has to send indices as global ordinals, if
    // the source and target graphs do not have the same column Map.
    // If LO doesn't match, the graphs could communicate using global
    // indices.  It could be possible that Node affects the graph's
    // storage format, but packAndPrepare should assume a common
    // communication format in any case.
    const row_graph_type* srcRowGraph = dynamic_cast<const row_graph_type*> (&source);
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      srcRowGraph == NULL, std::invalid_argument,
      ": The source object must be a RowGraph with matching first three "
      "template parameters.");

    // If the source object is actually a CrsGraph, we can use view
    // mode instead of copy mode to access the entries in each row,
    // if the graph is not fill complete.
    const this_type* srcCrsGraph = dynamic_cast<const this_type*> (&source);

    const map_type& srcRowMap = * (srcRowGraph->getRowMap ());
    const map_type& tgtRowMap = * (this->getRowMap ());
    const bool src_filled = srcRowGraph->isFillComplete ();
    Array<GO> row_copy;
    LO myid = 0;

    //
    // "Copy" part of "copy and permute."
    //
    if (src_filled || srcCrsGraph == NULL) {
      // If the source graph is fill complete, we can't use view mode,
      // because the data might be stored in a different format not
      // compatible with the expectations of view mode.  Also, if the
      // source graph is not a CrsGraph, we can't use view mode,
      // because RowGraph only provides copy mode access to the data.
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        size_t row_length = srcRowGraph->getNumEntriesInGlobalRow (gid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph->getGlobalRowCopy (gid, row_copy (), check_row_length);
        this->insertGlobalIndices (gid, row_copy ());
      }
    } else {
      for (size_t i = 0; i < numSameIDs; ++i, ++myid) {
        const GO gid = srcRowMap.getGlobalElement (myid);
        ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (gid, row);
        this->insertGlobalIndices (gid, row);
      }
    }

    //
    // "Permute" part of "copy and permute."
    //
    if (src_filled || srcCrsGraph == NULL) {
      for (LO i = 0; i < permuteToLIDs.size (); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs[i]);
        size_t row_length = srcRowGraph->getNumEntriesInGlobalRow (srcgid);
        row_copy.resize (row_length);
        size_t check_row_length = 0;
        srcRowGraph->getGlobalRowCopy (srcgid, row_copy (), check_row_length);
        this->insertGlobalIndices (mygid, row_copy ());
      }
    } else {
      for (LO i = 0; i < permuteToLIDs.size (); ++i) {
        const GO mygid = tgtRowMap.getGlobalElement (permuteToLIDs[i]);
        const GO srcgid = srcRowMap.getGlobalElement (permuteFromLIDs[i]);
        ArrayView<const GO> row;
        srcCrsGraph->getGlobalRowView (srcgid, row);
        this->insertGlobalIndices (mygid, row);
      }
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  packAndPrepare (const SrcDistObject& source,
                  const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                  Teuchos::Array<GlobalOrdinal> &exports,
                  const Teuchos::ArrayView<size_t> & numPacketsPerLID,
                  size_t& constantNumPackets,
                  Distributor &distor)
  {
    using Teuchos::Array;
    const char tfecfFuncName[] = "packAndPrepare";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": exportLIDs and numPacketsPerLID must have the same size.");
    typedef RowGraph<LocalOrdinal, GlobalOrdinal, Node> row_graph_type;
    const row_graph_type& srcGraph = dynamic_cast<const row_graph_type&> (source);

    // We don't check whether src_graph has had fillComplete called,
    // because it doesn't matter whether the *source* graph has been
    // fillComplete'd. The target graph can not be fillComplete'd yet.
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      this->isFillComplete (), std::runtime_error,
      ": The target graph of an Import or Export must not be fill complete.");

    srcGraph.pack (exportLIDs, exports, numPacketsPerLID, constantNumPackets, distor);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  pack (const Teuchos::ArrayView<const LocalOrdinal>& exportLIDs,
        Teuchos::Array<GlobalOrdinal>& exports,
        const Teuchos::ArrayView<size_t>& numPacketsPerLID,
        size_t& constantNumPackets,
        Distributor& distor) const
  {
    using Teuchos::Array;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    const char tfecfFuncName[] = "pack";
    (void) distor; // forestall "unused argument" compiler warning

    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      exportLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": exportLIDs and numPacketsPerLID must have the same size.");

    const map_type& srcMap = * (this->getMap ());
    constantNumPackets = 0;

    // Set numPacketsPerLID[i] to the number of entries owned by the
    // calling process in (local) row exportLIDs[i] of the graph, that
    // the caller wants us to send out.  Compute the total number of
    // packets (that is, entries) owned by this process in all the
    // rows that the caller wants us to send out.
    size_t totalNumPackets = 0;
    Array<GO> row;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      numPacketsPerLID[i] = row_length;
      totalNumPackets += row_length;
    }

    exports.resize (totalNumPackets);

    // Loop again over the rows to export, and pack rows of indices
    // into the output buffer.
    size_t exportsOffset = 0;
    for (LO i = 0; i < exportLIDs.size (); ++i) {
      const GO GID = srcMap.getGlobalElement (exportLIDs[i]);
      size_t row_length = this->getNumEntriesInGlobalRow (GID);
      row.resize (row_length);
      size_t check_row_length = 0;
      this->getGlobalRowCopy (GID, row (), check_row_length);
      typename Array<GO>::const_iterator row_iter = row.begin();
      typename Array<GO>::const_iterator row_end = row.end();
      size_t j = 0;
      for (; row_iter != row_end; ++row_iter, ++j) {
        exports[exportsOffset+j] = *row_iter;
      }
      exportsOffset += row.size();
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  unpackAndCombine (const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                    const Teuchos::ArrayView<const GlobalOrdinal> &imports,
                    const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                    size_t constantNumPackets,
                    Distributor& /* distor */,
                    CombineMode /* CM */)
  {
    // FIXME (mfh 02 Apr 2012) REPLACE combine mode has a perfectly
    // reasonable meaning, whether or not the matrix is fill complete.
    // It's just more work to implement.

    // We are not checking the value of the CombineMode input
    // argument.  For CrsGraph, we only support import/export
    // operations if fillComplete has not yet been called.  Any
    // incoming column-indices are inserted into the target graph. In
    // this context, CombineMode values of ADD vs INSERT are
    // equivalent. What is the meaning of REPLACE for CrsGraph? If a
    // duplicate column-index is inserted, it will be compressed out
    // when fillComplete is called.
    //
    // Note: I think REPLACE means that an existing row is replaced by
    // the imported row, i.e., the existing indices are cleared. CGB,
    // 6/17/2010

    const char tfecfFuncName[] = "unpackAndCombine";
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      importLIDs.size() != numPacketsPerLID.size(), std::runtime_error,
      ": importLIDs and numPacketsPerLID must have the same size.");
    TEUCHOS_TEST_FOR_EXCEPTION_CLASS_FUNC(
      isFillComplete (), std::runtime_error,
      ": Import or Export operations are not allowed on the destination "
      "CrsGraph if it is fill complete.");
    size_t importsOffset = 0;

    typedef typename ArrayView<const LocalOrdinal>::const_iterator iter_type;
    iter_type impLIDiter = importLIDs.begin();
    iter_type impLIDend = importLIDs.end();

    for (size_t i = 0; impLIDiter != impLIDend; ++impLIDiter, ++i) {
      LocalOrdinal row_length = numPacketsPerLID[i];
      ArrayView<const GlobalOrdinal> row (&imports[importsOffset], row_length);
      insertGlobalIndicesFiltered (this->getMap ()->getGlobalElement (*impLIDiter), row);
      importsOffset += row_length;
    }
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void
  CrsGraph<LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> ,  typename KokkosClassic::DefaultKernels<void,LocalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  removeEmptyProcessesInPlace (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& newMap)
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // We'll set all the state "transactionally," so that this method
    // satisfies the strong exception guarantee.  This object's state
    // won't be modified until the end of this method.
    RCP<const map_type> rowMap, domainMap, rangeMap, colMap;
    RCP<import_type> importer;
    RCP<export_type> exporter;

    rowMap = newMap;
    RCP<const Comm<int> > newComm =
      (newMap.is_null ()) ? null : newMap->getComm ();

    if (! domainMap_.is_null ()) {
      if (domainMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original domain and row Maps are identical.
        // In that case, we need only replace the original domain Map
        // with the new Map.  This ensures that the new domain and row
        // Maps _stay_ identical.
        domainMap = newMap;
      } else {
        domainMap = domainMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! rangeMap_.is_null ()) {
      if (rangeMap_.getRawPtr () == rowMap_.getRawPtr ()) {
        // Common case: original range and row Maps are identical.  In
        // that case, we need only replace the original range Map with
        // the new Map.  This ensures that the new range and row Maps
        // _stay_ identical.
        rangeMap = newMap;
      } else {
        rangeMap = rangeMap_->replaceCommWithSubset (newComm);
      }
    }
    if (! colMap.is_null ()) {
      colMap = colMap_->replaceCommWithSubset (newComm);
    }

    // (Re)create the Export and / or Import if necessary.
    if (! newComm.is_null ()) {
      RCP<ParameterList> params = this->getNonconstParameterList (); // could be null
      //
      // The operations below are collective on the new communicator.
      //
      // (Re)create the Export object if necessary.  If I haven't
      // called fillComplete yet, I don't have a rangeMap, so I must
      // first check if the _original_ rangeMap is not null.  Ditto
      // for the Import object and the domain Map.
      if (! rangeMap_.is_null () &&
          rangeMap != rowMap &&
          ! rangeMap->isSameAs (*rowMap)) {
        if (params.is_null () || ! params->isSublist ("Export")) {
          exporter = rcp (new export_type (rowMap, rangeMap));
        }
        else {
          RCP<ParameterList> exportSublist = sublist (params, "Export", true);
          exporter = rcp (new export_type (rowMap, rangeMap, exportSublist));
        }
      }
      // (Re)create the Import object if necessary.
      if (! domainMap_.is_null () &&
          domainMap != colMap &&
          ! domainMap->isSameAs (*colMap)) {
        if (params.is_null () || ! params->isSublist ("Import")) {
          importer = rcp (new import_type (domainMap, colMap));
        } else {
          RCP<ParameterList> importSublist = sublist (params, "Import", true);
          importer = rcp (new import_type (domainMap, colMap, importSublist));
        }
      }
    } // if newComm is not null

    // Defer side effects until the end.  If no destructors throw
    // exceptions (they shouldn't anyway), then this method satisfies
    // the strong exception guarantee.
    exporter_ = exporter;
    importer_ = importer;
    rowMap_ = rowMap;
    // mfh 31 Mar 2013: DistObject's map_ is the row Map of a CrsGraph
    // or CrsMatrix.  CrsGraph keeps a redundant pointer (rowMap_) to
    // the same object.  We might want to get rid of this redundant
    // pointer sometime, but for now, we'll leave it alone and just
    // set map_ to the same object.
    this->map_ = rowMap;
    domainMap_ = domainMap;
    rangeMap_ = rangeMap;
    colMap_ = colMap;
  }


  template <class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  bool
  CrsGraph<LocalOrdinal, GlobalOrdinal,
           Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>,
           typename KokkosClassic::DefaultKernels<void, LocalOrdinal,
             Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType> >::SparseOps>::
  hasRowInfo () const
  {
#ifdef HAVE_TPETRA_DEBUG
    bool actuallyHasRowInfo = true;
    if (indicesAreAllocated() && getProfileType() == StaticProfile && rowPtrs_ == null) {
      actuallyHasRowInfo = false;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(
      actuallyHasRowInfo != haveRowInfo_,
      std::logic_error, "Internal logic error. Please contact Tpetra team.");
#endif // HAVE_TPETRA_DEBUG
    return haveRowInfo_;
  }


} // namespace Tpetra




#endif // TPETRA_CRSGRAPH_DEF_HPP
