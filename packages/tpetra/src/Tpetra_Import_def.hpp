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

#ifndef TPETRA_IMPORT_DEF_HPP
#define TPETRA_IMPORT_DEF_HPP

#ifdef DOXYGEN_USE_ONLY
#  include <Tpetra_Import_decl.hpp>
#endif // DOXYGEN_USE_ONLY

#include <Tpetra_Distributor.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_ImportExportData.hpp>
#include <Tpetra_Util.hpp>
#include <Teuchos_as.hpp>

namespace Tpetra {
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source,
          const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target) :
    debug_ (false)
  {
    using Teuchos::rcp;
    using std::cerr;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    if (debug_) {
      cerr << source->getComm ()->getRank () << ": Import ctor" << endl;
    }
    ImportData_ = rcp (new data_type (source, target));
    setupSamePermuteRemote ();
    if (debug_) {
      cerr << source->getComm ()->getRank ()
           << ": Import ctor: setupSamePermuteRemote done" << endl;
    }
    if (source->isDistributed ()) {
      setupExport ();
    }
    if (debug_) {
      cerr << source->getComm ()->getRank ()
           << ": Import ctor: setupExport done" << endl;
    }
    remoteGIDs_ = null; // Don't need these anymore.
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & source,
          const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > & target,
          const Teuchos::RCP<Teuchos::ParameterList>& plist) :
    debug_ (false)
  {
    using Teuchos::rcp;
    using std::cerr;
    using std::endl;
    typedef ImportExportData<LocalOrdinal,GlobalOrdinal,Node> data_type;

    // Read "debug" parameter from the input ParameterList.
    bool debug = false;
    try {
      debug = plist->get<bool> ("debug");
    } catch (Teuchos::Exceptions::InvalidParameter&) {}
    debug_ = debug;

    if (debug_) {
      cerr << source->getComm ()->getRank () << ": Import ctor" << endl;
    }
    ImportData_ = rcp (new data_type (source, target, plist));
    setupSamePermuteRemote ();
    if (debug_) {
      cerr << source->getComm ()->getRank ()
           << ": Import ctor: setupSamePermuteRemote done" << endl;
    }
    if (source->isDistributed ()) {
      setupExport ();
    }
    if (debug_) {
      cerr << source->getComm ()->getRank ()
           << ": Import ctor: setupExport done" << endl;
    }
    remoteGIDs_ = null; // Don't need this anymore.
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  Import (const Import<LocalOrdinal,GlobalOrdinal,Node>& rhs)
    : ImportData_ (rhs.ImportData_)
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>::~Import()
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumSameIDs() const {
    return ImportData_->numSameIDs_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumPermuteIDs() const {
    return ImportData_->permuteFromLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteFromLIDs() const {
    return ImportData_->permuteFromLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getPermuteToLIDs() const {
    return ImportData_->permuteToLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumRemoteIDs() const {
    return ImportData_->remoteLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getRemoteLIDs() const {
    return ImportData_->remoteLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  size_t Import<LocalOrdinal,GlobalOrdinal,Node>::getNumExportIDs() const {
    return ImportData_->exportLIDs_.size();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const LocalOrdinal>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportLIDs() const {
    return ImportData_->exportLIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  ArrayView<const int>
  Import<LocalOrdinal,GlobalOrdinal,Node>::getExportImageIDs() const {
    return ImportData_->exportImageIDs_();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getSourceMap() const {
    return ImportData_->source_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  const RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getTargetMap() const {
    return ImportData_->target_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Distributor &
  Import<LocalOrdinal,GlobalOrdinal,Node>::getDistributor() const {
    return ImportData_->distributor_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Import<LocalOrdinal,GlobalOrdinal,Node>&
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  operator= (const Import<LocalOrdinal,GlobalOrdinal,Node>& rhs) {
    if (&rhs != this) {
      ImportData_ = rhs.ImportData_;
    }
    return *this;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::
  print (std::ostream& os) const
  {
    using Teuchos::Comm;
    using Teuchos::getFancyOStream;
    using Teuchos::RCP;
    using Teuchos::rcpFromRef;
    using Teuchos::toString;
    using std::endl;

    RCP<const Comm<int> > comm = getSourceMap ()->getComm ();
    const int myImageID = comm->getRank ();
    const int numImages = comm->getSize ();
    for (int imageCtr = 0; imageCtr < numImages; ++imageCtr) {
      if (myImageID == imageCtr) {
        os << endl;
        if (myImageID == 0) { // I'm the root node (only output this info once)
          os << "Import Data Members:" << endl;
        }
        os << "Image ID       : " << myImageID << endl;

        os << "permuteFromLIDs: " << toString (getPermuteFromLIDs ()) << endl;
        os << "permuteToLIDs  : " << toString (getPermuteToLIDs ()) << endl;
        os << "remoteLIDs     : " << toString (getRemoteLIDs ()) << endl;
        os << "exportLIDs     : " << toString (getExportLIDs ()) << endl;
        os << "exportImageIDs : " << toString (getExportImageIDs ()) << endl;

        os << "numSameIDs     : " << getNumSameIDs () << endl;
        os << "numPermuteIDs  : " << getNumPermuteIDs () << endl;
        os << "numRemoteIDs   : " << getNumRemoteIDs () << endl;
        os << "numExportIDs   : " << getNumExportIDs () << endl;
      }
      // A few global barriers give output a chance to complete.
      comm->barrier();
      comm->barrier();
      comm->barrier();
    }

    const bool printMaps = false;
    if (printMaps) {
      if (myImageID == 0) {
        os << endl << endl << "Source Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getSourceMap();
      comm->barrier();

      if (myImageID == 0) {
        os << endl << endl << "Target Map:" << endl << std::flush;
      }
      comm->barrier();
      os << *getTargetMap();
      comm->barrier();
    }

    // It's also helpful for debugging to print the Distributor
    // object.  Epetra_Import::Print() does this, so we can do a
    // side-by-side comparison.
    if (myImageID == 0) {
      os << endl << endl << "Distributor:" << endl << std::flush;
    }
    comm->barrier();
    getDistributor().describe (*(getFancyOStream (rcpFromRef (os))),
                               Teuchos::VERB_EXTREME);
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Import<LocalOrdinal,GlobalOrdinal,Node>::
  setupSamePermuteRemote()
  {
    const Map<LocalOrdinal,GlobalOrdinal,Node> & source = *getSourceMap();
    const Map<LocalOrdinal,GlobalOrdinal,Node> & target = *getTargetMap();
    ArrayView<const GlobalOrdinal> sourceGIDs = source.getNodeElementList();
    ArrayView<const GlobalOrdinal> targetGIDs = target.getNodeElementList();

    const bool fasterImplementation = false;
    if (fasterImplementation) {
      typedef typename ArrayView<const GlobalOrdinal>::size_type size_type;
      const LocalOrdinal LINVALID = Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();

      const GlobalOrdinal* const rawSrcGids = sourceGIDs.getRawPtr ();
      const GlobalOrdinal* const rawTgtGids = targetGIDs.getRawPtr ();
      const size_type numSrcGids = sourceGIDs.size ();
      const size_type numTgtGids = targetGIDs.size ();
      const size_type numGids = std::min (numSrcGids, numTgtGids);

      size_type numSameGids = 0;
      for ( ; numSameGids < numGids && rawSrcGids[numSameGids] == rawTgtGids[numSameGids]; ++numSameGids)
        {} // third clause of 'for' does everything
      ImportData_->numSameIDs_ = numSameGids;

      remoteGIDs_ = rcp (new Array<GlobalOrdinal> ());
      Array<GlobalOrdinal>& remoteGids = *remoteGIDs_;
      Array<LocalOrdinal>& permuteToLIDs = ImportData_->permuteToLIDs_;
      Array<LocalOrdinal>& permuteFromLIDs = ImportData_->permuteFromLIDs_;
      Array<LocalOrdinal>& remoteLIDs = ImportData_->remoteLIDs_;

      for (size_type targetIndex = numSameGids; targetIndex < numTgtGids; ++targetIndex) {
        const GlobalOrdinal curTargetGid = rawTgtGids[targetIndex];
        // getLocalElement() returns LINVALID if the GID isn't in the source Map.
        // This saves us a lookup (which isNodeGlobalElement() would do).
        const LocalOrdinal srcLid = source.getLocalElement (curTargetGid);
        const LocalOrdinal tgtLid = target.getLocalElement (curTargetGid);

        if (srcLid != LINVALID) { // if source.isNodeGlobalElement (curTargetGid)
          permuteToLIDs.push_back (tgtLid);
          permuteFromLIDs.push_back (srcLid);
        } else {
          remoteGids.push_back (curTargetGid);
          remoteLIDs.push_back (tgtLid);
        }
      }
    } else {
      // Compute numSameIDs_:
      //
      // Iterate through the source and target GID lists.  If the i-th
      // GID of both is the same, increment numSameIDs_ and try the
      // next.  As soon as you come to a nonmatching pair, give up.
      //
      // The point of numSameIDs_ is for the common case of an Import
      // where all the overlapping GIDs are at the end of the source
      // Map, but otherwise the source and target Maps are the same.
      // This allows a fast contiguous copy for the initial "same IDs."
      typename ArrayView<const GlobalOrdinal>::iterator sourceIter = sourceGIDs.begin(),
        targetIter = targetGIDs.begin();
      while (sourceIter != sourceGIDs.end() &&
             targetIter != targetGIDs.end() &&
             *sourceIter == *targetIter) {
        ++ImportData_->numSameIDs_;
        ++sourceIter;
        ++targetIter;
      }
      // targetIter should now point either to the GID of the first
      // non-same entry in targetGIDs, or to the end of targetGIDs (if
      // all the entries were the same).

      // Compute IDs to be permuted, vs. IDs to be received ("imported";
      // called "remote" IDs).
      //
      // IDs to permute are in both the source and target Maps, which
      // means we don't have to send or receive them, but we do have to
      // rearrange (permute) them in general.  (We've already identified
      // an initial stretch of IDs which can be copied without
      // rearrangement; the iterator targetIter is past that point.)
      // IDs to receive are in the target Map, but not the source Map.
      //
      // How do the following code and its equivalent in Export differ?
      //
      // 1. Export uses sourceIter, whereas Import uses targetIter.
      //
      // 2. Import collects remoteGIDs_ (target Map GIDs that are not in
      //    the source Map), which is a separate array.  Export can use
      //    exportGIDs_, which is an array belonging to its
      //    ImportExportData object.
      remoteGIDs_ = rcp( new Array<GlobalOrdinal>() );
      for (; targetIter != targetGIDs.end(); ++targetIter) {
        const GlobalOrdinal curTargetGID = *targetIter;
        if (source.isNodeGlobalElement (curTargetGID)) {
          // The current process owns this GID, for both the source and
          // the target Maps.  Determine the LIDs for this GID on both
          // Maps and add them to the permutation lists.
          ImportData_->permuteToLIDs_.push_back (target.getLocalElement (curTargetGID));
          ImportData_->permuteFromLIDs_.push_back (source.getLocalElement (curTargetGID));
        }
        else {
          // The current GID is owned by this process in the target Map,
          // but is not owned by this process in the source Map.  That
          // means the Import operation has to receive it from another
          // process.  Store it in the "remote" (incoming) list, along
          // with its destination LID on this process.
          //
          // remoteLIDs_ is the list of this process' LIDs that it has
          // to receive from other processes.  Since this is an Import,
          // and therefore the source Map is nonoverlapping, we know
          // that each remote LID can receive from only one process.
          remoteGIDs_->push_back (curTargetGID);
          ImportData_->remoteLIDs_.push_back (target.getLocalElement (curTargetGID));
        }
      }
    } // if using the hopefully faster implementation

    TPETRA_ABUSE_WARNING(
      getNumRemoteIDs() > 0 && ! source.isDistributed(),
      std::runtime_error,
      "::setupSamePermuteRemote(): Target has remote LIDs but Source is not "
      "distributed globally." << std::endl
      << "Importing to a submap of the target map.");
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void Import<LocalOrdinal,GlobalOrdinal,Node>::setupExport() {
    using Teuchos::arcp;
    using Teuchos::Array;
    using Teuchos::ArrayRCP;
    using Teuchos::ArrayView;
    using Teuchos::null;
    using std::cerr;
    using std::endl;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef typename Array<int>::difference_type size_type;
    const Map<LO, GO, Node> & source = *getSourceMap ();

    if (debug_) {
      const int myRank = source.getComm ()->getRank ();
      cerr << myRank << ": Import::setupExport" << endl;
    }

    // For each entry remoteGIDs[i], remoteProcIDs[i] will contain
    // the process ID of the process that owns that GID.
    ArrayView<GO> remoteGIDs = (*remoteGIDs_) ();
    Array<int> remoteProcIDs (remoteGIDs.size ());

    // lookup == IDNotPresent means that the source Map wasn't able to
    // figure out to which processes one or more of the GIDs in the
    // given list of remoteGIDs belong.
    //
    // The previous abuse warning said "The target Map has GIDs not
    // found in the source Map."  This statement could be confusing,
    // because it doesn't refer to ownership by the current process,
    // but rather to ownership by _any_ process participating in the
    // Map.  (It could not possibly refer to ownership by the current
    // process, since remoteGIDs is exactly the list of GIDs owned by
    // the target Map but not owned by the source Map.  It was
    // constructed that way by setupSamePermuteRemote().)
    //
    // What this statement means is that the source and target Maps
    // don't contain the same set of GIDs globally (over all
    // processes).  That is, there is at least one GID owned by some
    // process in the target Map, which is not owned by _any_ process
    // in the source Map.
    const LookupStatus lookup =
      source.getRemoteIndexList (remoteGIDs, remoteProcIDs ());
    if (debug_) {
      const int myRank = source.getComm ()->getRank ();
      cerr << myRank << ": Import::setupExport: finished lookup" << endl;
    }

    TPETRA_ABUSE_WARNING( lookup == IDNotPresent, std::runtime_error,
      "::setupExport(): the source Map wasn't able to figure out which process "
      "owns one or more of the GIDs in the list of remote GIDs.  This probably "
      "means that there is at least one GID owned by some process in the target"
      " Map which is not owned by any process in the source Map.  (That is, the"
      " source and target Maps do not contain the same set of GIDs globally.)");

    // Ignore remote GIDs that aren't owned by any process in the
    // source Map.  getRemoteIndexList() gives each of these a process
    // ID of -1.
    if (lookup == IDNotPresent) {
      const size_type numInvalidRemote =
        std::count_if (remoteProcIDs.begin (), remoteProcIDs.end (),
                       std::bind1st (std::equal_to<int> (), -1));
      // If all of them are invalid, we can delete the whole array.
      const size_type totalNumRemote = getNumRemoteIDs ();
      if (numInvalidRemote == totalNumRemote) {
        // all remotes are invalid; we have no remotes; we can delete the remotes
        remoteProcIDs.clear ();
        (*remoteGIDs_).clear (); // This invalidates the view remoteGIDs
        ImportData_->remoteLIDs_.clear();
      }
      else {
        // Some remotes are valid; we need to keep the valid ones.
        // Pack and resize remoteProcIDs, remoteGIDs_, and
        // remoteLIDs_.
        size_type numValidRemote = 0;
#ifdef HAVE_TPETRA_DEBUG
        ArrayView<GlobalOrdinal> remoteGIDsPtr = remoteGIDs;
#else
        GlobalOrdinal* const remoteGIDsPtr = remoteGIDs.getRawPtr ();
#endif // HAVE_TPETRA_DEBUG
        for (size_type r = 0; r < totalNumRemote; ++r) {
          // Pack in all the valid remote PIDs and GIDs.
          if (remoteProcIDs[r] != -1) {
            remoteProcIDs[numValidRemote] = remoteProcIDs[r];
            remoteGIDsPtr[numValidRemote] = remoteGIDsPtr[r];
            ImportData_->remoteLIDs_[numValidRemote] = ImportData_->remoteLIDs_[r];
            ++numValidRemote;
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          numValidRemote != totalNumRemote - numInvalidRemote, std::logic_error,
          "Tpetra::Import::setupExport(): After removing invalid remote GIDs and"
          " packing the valid remote GIDs, numValidRemote = " << numValidRemote
          << " != totalNumRemote - numInvalidRemote = "
          << totalNumRemote - numInvalidRemote
          << ".  Please report this bug to the Tpetra developers.");

        remoteProcIDs.resize (numValidRemote);
        (*remoteGIDs_).resize (numValidRemote);
        ImportData_->remoteLIDs_.resize (numValidRemote);
      }
      // Revalidate the view after clear or resize.
      remoteGIDs = (*remoteGIDs_)();
    }

    // Sort remoteProcIDs in ascending order, and apply the resulting
    // permutation to remoteGIDs_ and remoteLIDs_.  This ensures that
    // remoteProcIDs[i], remoteGIDs_[i], and remoteLIDs_[i] all refer
    // to the same thing.
    sort3 (remoteProcIDs.begin(),
           remoteProcIDs.end(),
           remoteGIDs.begin(),
           ImportData_->remoteLIDs_.begin());

    // Call the Distributor's createFromRecvs() method to turn the
    // remote GIDs and their owning processes into a send-and-receive
    // communication plan.  remoteGIDs and remoteProcIDs_ are input;
    // exportGIDs and exportProcIDs_ are output arrays which are
    // allocated by createFromRecvs().
    ArrayRCP<GO> exportGIDs;
    ImportData_->distributor_.createFromRecvs (remoteGIDs ().getConst (),
                                               remoteProcIDs, exportGIDs,
                                               ImportData_->exportImageIDs_);
    if (debug_) {
      const int myRank = source.getComm ()->getRank ();
      cerr << myRank << ": Import::setupExport: finished createFromRecvs" << endl;
    }

    // Find the LIDs corresponding to the (outgoing) GIDs in
    // exportGIDs.  For sparse matrix-vector multiply, this tells the
    // calling process how to index into the source vector to get the
    // elements which it needs to send.
    const size_type numExportIDs = exportGIDs.size ();
    if (numExportIDs > 0) {
      ImportData_->exportLIDs_ = arcp<LO> (numExportIDs);

      ArrayView<const GO> expGIDs = exportGIDs ();
      ArrayView<LO> expLIDs = ImportData_->exportLIDs_ ();
      for (size_type k = 0; k < numExportIDs; ++k) {
        expLIDs[k] = source.getLocalElement (expGIDs[k]);
      }
    }
  }
} // namespace Tpetra

// Explicit instantiation macro.
// Only invoke this when in the Tpetra namespace.
// Most users do not need to use this.
//
// LO: The local ordinal type.
// GO: The global ordinal type.
// NODE: The Kokkos Node type.
#define TPETRA_IMPORT_INSTANT(LO, GO, NODE) \
  \
  template class Import< LO , GO , NODE >;

#endif // TPETRA_IMPORT_DEF_HPP
