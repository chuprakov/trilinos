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

/// \file Tpetra_Map_def.hpp
///
/// Implementation of the methods of Tpetra::Map, and of related
/// nonmember constructors for Tpetra::Map.

#ifndef TPETRA_MAP_DEF_HPP
#define TPETRA_MAP_DEF_HPP

#include <Tpetra_Directory.hpp> // must include for implicit instantiation to work
#include <Tpetra_Details_FixedHashTable.hpp>
#include <Tpetra_Util.hpp>
#include <Teuchos_as.hpp>
#include <stdexcept>

#ifdef DOXYGEN_USE_ONLY
  #include "Tpetra_Map_decl.hpp"
#endif

namespace Tpetra {
  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map () :
    indexBase_ (0),
    numGlobalElements_ (0),
    numLocalElements_ (0),
    minMyGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    maxMyGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    minAllGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    maxAllGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    firstContiguousGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    lastContiguousGID_ (Teuchos::OrdinalTraits<GlobalOrdinal>::invalid ()),
    uniform_ (false), // trivially
    contiguous_ (false),
    distributed_ (false) // no communicator yet
  {}

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements,
       GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
       LocalGlobal lOrG,
       const Teuchos::RCP<Node> &node) :
    comm_ (comm),
    node_ (node),
    uniform_ (true)
  {
    using Teuchos::as;
    using Teuchos::broadcast;
    using Teuchos::outArg;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_MAX;
    using Teuchos::typeName;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;
    const GST GSTI = Teuchos::OrdinalTraits<GST>::invalid ();

#ifdef HAVE_TPETRA_DEBUG
    // In debug mode only, check whether numGlobalElements and
    // indexBase are the same over all processes in the communicator.
    {
      GST proc0NumGlobalElements = numGlobalElements;
      broadcast<int, GST> (*comm_, 0, outArg (proc0NumGlobalElements));
      GST minNumGlobalElements = numGlobalElements;
      GST maxNumGlobalElements = numGlobalElements;
      reduceAll<int, GST> (*comm, REDUCE_MIN, numGlobalElements, outArg (minNumGlobalElements));
      reduceAll<int, GST> (*comm, REDUCE_MAX, numGlobalElements, outArg (maxNumGlobalElements));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minNumGlobalElements != maxNumGlobalElements || numGlobalElements != minNumGlobalElements,
        std::invalid_argument,
        "Tpetra::Map constructor: All processes must provide the same number "
        "of global elements.  Process 0 set numGlobalElements = "
        << proc0NumGlobalElements << ".  The calling process "
        << comm->getRank () << " set numGlobalElements = " << numGlobalElements
        << ".  The min and max values over all processes are "
        << minNumGlobalElements << " resp. " << maxNumGlobalElements << ".");

      GO proc0IndexBase = indexBase;
      broadcast<int, GO> (*comm_, 0, outArg (proc0IndexBase));
      GO minIndexBase = indexBase;
      GO maxIndexBase = indexBase;
      reduceAll<int, GO> (*comm, REDUCE_MIN, indexBase, outArg (minIndexBase));
      reduceAll<int, GO> (*comm, REDUCE_MAX, indexBase, outArg (maxIndexBase));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minIndexBase != maxIndexBase || indexBase != minIndexBase,
        std::invalid_argument,
        "Tpetra::Map constructor: "
        "All processes must provide the same indexBase argument.  "
        "Process 0 set indexBase = " << proc0IndexBase << ".  The calling "
        "process " << comm->getRank () << " set indexBase = " << indexBase
        << ".  The min and max values over all processes are "
        << minIndexBase << " resp. " << maxIndexBase << ".");
    }
#endif // HAVE_TPETRA_DEBUG

    // Distribute the elements across the processes in the given
    // communicator so that global IDs (GIDs) are
    //
    // - Nonoverlapping (only one process owns each GID)
    // - Contiguous (the sequence of GIDs is nondecreasing, and no two
    //   adjacent GIDs differ by more than one)
    // - As evenly distributed as possible (the numbers of GIDs on two
    //   different processes do not differ by more than one)

    // All processes have the same numGlobalElements, but we still
    // need to check that it is valid.  numGlobalElements must be
    // positive and not the "invalid" value (GSTI).
    //
    // This comparison looks funny, but it avoids compiler warnings
    // for comparing unsigned integers (numGlobalElements_in is a
    // GST, which is unsigned) while still working if we
    // later decide to make GST signed.
    TEUCHOS_TEST_FOR_EXCEPTION(
      (numGlobalElements < 1 && numGlobalElements != 0),
      std::invalid_argument,
      "Tpetra::Map constructor: numGlobalElements (= "
      << numGlobalElements << ") must be nonnegative.");

    TEUCHOS_TEST_FOR_EXCEPTION(
      numGlobalElements == GSTI, std::invalid_argument,
      "Tpetra::Map constructor: You provided numGlobalElements = Teuchos::"
      "OrdinalTraits<Tpetra::global_size_t>::invalid().  This version of the "
      "constructor requires a valid value of numGlobalElements.  You "
      "probably mistook this constructor for the \"contiguous nonuniform\" "
      "constructor, which can compute the global number of elements for you "
      "if you set numGlobalElements to that value.");

    size_t numLocalElements = 0; // will set below
    if (lOrG == GloballyDistributed) {
      // Compute numLocalElements:
      //
      // If numGlobalElements == numProcs * B + remainder,
      // then Proc r gets B+1 elements if r < remainder,
      // and B elements if r >= remainder.
      //
      // This strategy is valid for any value of numGlobalElements and
      // numProcs, including the following border cases:
      //   - numProcs == 1
      //   - numLocalElements < numProcs
      //
      // In the former case, remainder == 0 && numGlobalElements ==
      // numLocalElements.  In the latter case, remainder ==
      // numGlobalElements && numLocalElements is either 0 or 1.
      const GST numProcs = as<GST> (comm_->getSize ());
      const GST myRank = as<GST> (comm_->getRank ());
      const GST quotient  = numGlobalElements / numProcs;
      const GST remainder = numGlobalElements - quotient * numProcs;

      GO startIndex;
      if (myRank < remainder) {
        numLocalElements = as<size_t> (1) + as<size_t> (quotient);
        // myRank was originally an int, so it should never overflow
        // reasonable GO types.
        startIndex = as<GO> (myRank) * as<GO> (numLocalElements);
      } else {
        numLocalElements = as<size_t> (quotient);
        startIndex = as<GO> (myRank) * as<GO> (numLocalElements) +
          as<GO> (remainder);
      }

      minMyGID_  = indexBase + startIndex;
      maxMyGID_  = indexBase + startIndex + numLocalElements - 1;
      minAllGID_ = indexBase;
      maxAllGID_ = indexBase + numGlobalElements - 1;
      distributed_ = (numProcs > 1);
    }
    else {  // lOrG == LocallyReplicated
      numLocalElements = as<size_t> (numGlobalElements);
      minMyGID_ = indexBase;
      maxMyGID_ = indexBase + numGlobalElements - 1;
      distributed_ = false;
    }

    minAllGID_ = indexBase;
    maxAllGID_ = indexBase + numGlobalElements - 1;
    indexBase_ = indexBase;
    numGlobalElements_ = numGlobalElements;
    numLocalElements_ = numLocalElements;
    firstContiguousGID_ = minMyGID_;
    lastContiguousGID_ = maxMyGID_;
    contiguous_ = true;

    // Create the Directory on demand in getRemoteIndexList().
    //setupDirectory ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements,
       size_t numLocalElements,
       GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
       const Teuchos::RCP<Node> &node) :
    comm_ (comm),
    node_ (node),
    uniform_ (false)
  {
    using Teuchos::as;
    using Teuchos::broadcast;
    using Teuchos::outArg;
    using Teuchos::reduceAll;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_SUM;
    using Teuchos::scan;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;
    const GST GSTI = Teuchos::OrdinalTraits<GST>::invalid ();

#ifdef HAVE_TPETRA_DEBUG
    // Keep this for later debug checks.
    GST debugGlobalSum = 0; // Will be global sum of numLocalElements
    reduceAll<int, GST> (*comm, REDUCE_SUM, as<GST> (numLocalElements),
                         outArg (debugGlobalSum));
    // In debug mode only, check whether numGlobalElements and
    // indexBase are the same over all processes in the communicator.
    {
      GST proc0NumGlobalElements = numGlobalElements;
      broadcast<int, GST> (*comm_, 0, outArg (proc0NumGlobalElements));
      GST minNumGlobalElements = numGlobalElements;
      GST maxNumGlobalElements = numGlobalElements;
      reduceAll<int, GST> (*comm, REDUCE_MIN, numGlobalElements, outArg (minNumGlobalElements));
      reduceAll<int, GST> (*comm, REDUCE_MAX, numGlobalElements, outArg (maxNumGlobalElements));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minNumGlobalElements != maxNumGlobalElements || numGlobalElements != minNumGlobalElements,
        std::invalid_argument,
        "Tpetra::Map constructor: All processes must provide the same number "
        "of global elements.  This is true even if that argument is Teuchos::"
        "OrdinalTraits<global_size_t>::invalid() to signal that the Map should "
        "compute the global number of elements.  Process 0 set numGlobalElements"
        " = " << proc0NumGlobalElements << ".  The calling process "
        << comm->getRank () << " set numGlobalElements = " << numGlobalElements
        << ".  The min and max values over all processes are "
        << minNumGlobalElements << " resp. " << maxNumGlobalElements << ".");

      GO proc0IndexBase = indexBase;
      broadcast<int, GO> (*comm_, 0, outArg (proc0IndexBase));
      GO minIndexBase = indexBase;
      GO maxIndexBase = indexBase;
      reduceAll<int, GO> (*comm, REDUCE_MIN, indexBase, outArg (minIndexBase));
      reduceAll<int, GO> (*comm, REDUCE_MAX, indexBase, outArg (maxIndexBase));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minIndexBase != maxIndexBase || indexBase != minIndexBase,
        std::invalid_argument,
        "Tpetra::Map constructor: "
        "All processes must provide the same indexBase argument.  "
        "Process 0 set indexBase = " << proc0IndexBase << ".  The calling "
        "process " << comm->getRank () << " set indexBase = " << indexBase
        << ".  The min and max values over all processes are "
        << minIndexBase << " resp. " << maxIndexBase << ".");

      // Make sure that the sum of numLocalElements over all processes
      // equals numGlobalElements.
      TEUCHOS_TEST_FOR_EXCEPTION(
        numGlobalElements != GSTI && debugGlobalSum != numGlobalElements,
        std::invalid_argument,
        "Tpetra::Map constructor: The sum of numLocalElements over all "
        "processes = " << debugGlobalSum << " != numGlobalElements = "
        << numGlobalElements << ".  If you would like this constructor to "
        "compute numGlobalElements for you, you may set numGlobalElements = "
        "Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid() on input.");
    }
#endif // HAVE_TPETRA_DEBUG

    // Distribute the elements across the nodes so that they are
    // - non-overlapping
    // - contiguous

    // This differs from the first Map constructor (that only takes a
    // global number of elements) in that the user has specified the
    // number of local elements, so that the elements are not
    // (necessarily) evenly distributed over the processes.

    // Compute my local offset.  This is an inclusive scan, so to get
    // the final offset, we subtract off the input.
    GO scanResult = 0;
    scan<int, GO> (*comm, REDUCE_SUM, numLocalElements, outArg (scanResult));
    const GO myOffset = scanResult - numLocalElements;

    if (numGlobalElements != GSTI) {
      numGlobalElements_ = numGlobalElements; // Use the user's value.
    } else {
      // Inclusive scan means that the last process has the final sum.
      // Rather than doing a reduceAll to get the sum of
      // numLocalElements, we can just have the last process broadcast
      // its result.  That saves us a round of log(numProcs) messages.
      const int numProcs = comm->getSize ();
      GST globalSum = scanResult;
      if (numProcs > 1) {
        broadcast (*comm, numProcs - 1, outArg (globalSum));
      }
      numGlobalElements_ = globalSum;

#ifdef HAVE_TPETRA_DEBUG
      // No need for an all-reduce here; both come from collectives.
      TEUCHOS_TEST_FOR_EXCEPTION(
        globalSum != debugGlobalSum, std::logic_error,
        "Tpetra::Map constructor (contiguous nonuniform): "
        "globalSum = " << globalSum << " != debugGlobalSum = " << debugGlobalSum
        << ".  Please report this bug to the Tpetra developers.");
#endif // HAVE_TPETRA_DEBUG
    }
    numLocalElements_ = numLocalElements;
    indexBase_ = indexBase;
    minAllGID_ = indexBase;
    // numGlobalElements might be GSTI; use numGlobalElements_;
    maxAllGID_ = indexBase + numGlobalElements_ - 1;
    minMyGID_ = indexBase + myOffset;
    maxMyGID_ = indexBase + myOffset + numLocalElements - 1;
    firstContiguousGID_ = minMyGID_;
    lastContiguousGID_ = maxMyGID_;
    contiguous_ = true;
    distributed_ = checkIsDist ();

    // Create the Directory on demand in getRemoteIndexList().
    //setupDirectory ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  Map (global_size_t numGlobalElements,
       const Teuchos::ArrayView<const GlobalOrdinal> &entryList,
       GlobalOrdinal indexBase,
       const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
       const Teuchos::RCP<Node> &node) :
    comm_ (comm),
    node_ (node),
    uniform_ (false)
  {
    using Teuchos::arcp;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::broadcast;
    using Teuchos::outArg;
    using Teuchos::ptr;
    using Teuchos::REDUCE_MAX;
    using Teuchos::REDUCE_MIN;
    using Teuchos::REDUCE_SUM;
    using Teuchos::reduceAll;
    using Teuchos::typeName;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef global_size_t GST;
    typedef typename ArrayView<const GO>::size_type size_type;
    const GST GSTI = Teuchos::OrdinalTraits<GST>::invalid ();

    // The user has specified the distribution of elements over the
    // processes, via entryList.  The distribution is not necessarily
    // contiguous or equally shared over the processes.

    // The length of entryList on this node is the number of local
    // elements (on this node), even though entryList contains global
    // indices.  We assume that the number of local elements can be
    // stored in a size_t; numLocalElements_ is a size_t, so this
    // variable and that should have the same type.
    const size_t numLocalElements = as<size_t> (entryList.size ());

#ifdef HAVE_TPETRA_DEBUG
    // Keep this for later debug checks.
    GST debugGlobalSum = 0; // Will be global sum of numLocalElements
    reduceAll<int, GST> (*comm, REDUCE_SUM, as<GST> (numLocalElements),
                         outArg (debugGlobalSum));
    // In debug mode only, check whether numGlobalElements and
    // indexBase are the same over all processes in the communicator.
    {
      GST proc0NumGlobalElements = numGlobalElements;
      broadcast<int, GST> (*comm_, 0, outArg (proc0NumGlobalElements));
      GST minNumGlobalElements = numGlobalElements;
      GST maxNumGlobalElements = numGlobalElements;
      reduceAll<int, GST> (*comm, REDUCE_MIN, numGlobalElements, outArg (minNumGlobalElements));
      reduceAll<int, GST> (*comm, REDUCE_MAX, numGlobalElements, outArg (maxNumGlobalElements));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minNumGlobalElements != maxNumGlobalElements || numGlobalElements != minNumGlobalElements,
        std::invalid_argument,
        "Tpetra::Map constructor: All processes must provide the same number "
        "of global elements.  This is true even if that argument is Teuchos::"
        "OrdinalTraits<global_size_t>::invalid() to signal that the Map should "
        "compute the global number of elements.  Process 0 set numGlobalElements"
        " = " << proc0NumGlobalElements << ".  The calling process "
        << comm->getRank () << " set numGlobalElements = " << numGlobalElements
        << ".  The min and max values over all processes are "
        << minNumGlobalElements << " resp. " << maxNumGlobalElements << ".");

      GO proc0IndexBase = indexBase;
      broadcast<int, GO> (*comm_, 0, outArg (proc0IndexBase));
      GO minIndexBase = indexBase;
      GO maxIndexBase = indexBase;
      reduceAll<int, GO> (*comm, REDUCE_MIN, indexBase, outArg (minIndexBase));
      reduceAll<int, GO> (*comm, REDUCE_MAX, indexBase, outArg (maxIndexBase));
      TEUCHOS_TEST_FOR_EXCEPTION(
        minIndexBase != maxIndexBase || indexBase != minIndexBase,
        std::invalid_argument,
        "Tpetra::Map constructor: "
        "All processes must provide the same indexBase argument.  "
        "Process 0 set indexBase = " << proc0IndexBase << ".  The calling "
        "process " << comm->getRank () << " set indexBase = " << indexBase
        << ".  The min and max values over all processes are "
        << minIndexBase << " resp. " << maxIndexBase << ".");

      // Make sure that the sum of numLocalElements over all processes
      // equals numGlobalElements.
      TEUCHOS_TEST_FOR_EXCEPTION(
        ((numGlobalElements != GSTI) && (debugGlobalSum != numGlobalElements)),
        std::invalid_argument,
        "Tpetra::Map constructor: The sum of entryList.size() over all "
        "processes = " << debugGlobalSum << " != numGlobalElements = "
        << numGlobalElements << ".  If you would like this constructor to "
        "compute numGlobalElements for you, you may set numGlobalElements = "
        "Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid() on input.");
    }
#endif // HAVE_TPETRA_DEBUG

    // FIXME (mfh 20 Feb 2013) The global reduction is redundant,
    // since the directory Map will have to do the same thing.  We
    // should actually do the scan and broadcast for the directory Map
    // here, and give the computed offsets to the directory Map's
    // constructor.
    if (numGlobalElements != GSTI) {
      numGlobalElements_ = numGlobalElements; // Use the user's value.
    } else { // The user wants us to compute the sum.
      reduceAll<int, GST> (*comm, REDUCE_SUM, as<GST> (numLocalElements),
                           outArg (numGlobalElements_));
    }

    // mfh 20 Feb 2013: We've never quite done the right thing for
    // duplicate GIDs here.  Duplicate GIDs have always been counted
    // distinctly in numLocalElements_, and thus should get a
    // different LID.  However, we've always used std::map or a hash
    // table for the GID -> LID lookup table, so distinct GIDs always
    // map to the same LID.  Furthermore, the order of the input GID
    // list matters, so it's not desirable to sort for determining
    // uniqueness.
    //
    // I've chosen for now to write this code as if the input GID list
    // contains no duplicates.  If this is not desired, we could use
    // the lookup table itself to determine uniqueness: If we haven't
    // seen the GID before, it gets a new LID and it's added to the
    // LID -> GID and GID -> LID tables.  If we have seen the GID
    // before, it doesn't get added to either table.  I would
    // implement this, but it would cost more to do the double lookups
    // in the table (one to check, and one to insert).

    numLocalElements_ = numLocalElements;
    indexBase_ = indexBase;

    minMyGID_ = indexBase_;
    maxMyGID_ = indexBase_;

    if (numLocalElements_ > 0) {
      // Find contiguous GID range, with the restriction that the
      // beginning of the range starts with the first entry.  While
      // doing so, fill in the LID -> GID table.
      lgMap_ = arcp<GO> (numLocalElements_);
      firstContiguousGID_ = entryList[0];
      lastContiguousGID_ = firstContiguousGID_+1;
      lgMap_[0] = firstContiguousGID_;
      size_t i = 1;
      for ( ; i < numLocalElements_; ++i) {
        const GO curGid = entryList[i];
        const LO curLid = as<LO> (i);

        if (lastContiguousGID_ != curGid) break;

        // Add the entry to the LID->GID table only after we know that
        // the current GID is in the initial contiguous sequence, so
        // that we don't repeat adding it in the first iteration of
        // the loop below over the remaining noncontiguous GIDs.
        lgMap_[curLid] = curGid;
        ++lastContiguousGID_;
      }
      --lastContiguousGID_;

      // [firstContiguousGID_, lastContigousGID_] is the initial
      // sequence of contiguous GIDs.  We can start the min and max
      // GID using this range.
      minMyGID_ = firstContiguousGID_;
      maxMyGID_ = lastContiguousGID_;

      // Compute the GID -> LID lookup table, _not_ including the
      // initial sequence of contiguous GIDs.
      ArrayView<const GO> nonContigEntries =
        entryList (as<size_type> (i), entryList.size () - as<size_type> (i));
      glMap_ = rcp (new global_to_local_table_type (nonContigEntries, as<LO> (i)));

      for ( ; i < numLocalElements_; ++i) {
        const GO curGid = entryList[i];
        const LO curLid = as<LO> (i);
        lgMap_[curLid] = curGid; // LID -> GID table

        // While iterating through entryList, we compute its
        // (process-local) min and max elements.
        if (curGid < minMyGID_) {
          minMyGID_ = curGid;
        }
        if (curGid > maxMyGID_) {
          maxMyGID_ = curGid;
        }
      }
    }
    else {
      // This insures tests for GIDs in the range
      // [firstContiguousGID_, lastContiguousGID_] fail for processes
      // with no local elements.
      firstContiguousGID_ = indexBase_+1;
      lastContiguousGID_ = indexBase_;
      glMap_ = rcp (new global_to_local_table_type (entryList)); // is empty
    }

    // Compute the min and max of all processes' GIDs.  If
    // numLocalElements_ == 0 on this process, minMyGID_ and maxMyGID_
    // are both indexBase_.  This is wrong, but fixing it would
    // require either a fancy sparse all-reduce, or a custom reduction
    // operator that ignores invalid values ("invalid" means
    // Teuchos::OrdinalTraits<GO>::invalid()).
    //
    // Also, while we're at it, use the same all-reduce to figure out
    // if the Map is distributed.  "Distributed" means that there is
    // at least one process with a number of local elements less than
    // the number of global elements.
    //
    // We're computing the min and max of all processes' GIDs using a
    // single MAX all-reduce, because min(x,y) = -max(-x,-y) (when x
    // and y are signed).  (This lets us combine the min and max into
    // a single all-reduce.)  If each process sets localDist=1 if its
    // number of local elements is strictly less than the number of
    // global elements, and localDist=0 otherwise, then a MAX
    // all-reduce on localDist tells us if the Map is distributed (1
    // if yes, 0 if no).  Thus, we can append localDist onto the end
    // of the data and get the global result from the all-reduce.
    if (std::numeric_limits<GO>::is_signed) {
      // Does my process NOT own all the elements?
      const GO localDist =
        (as<GST> (numLocalElements_) < numGlobalElements_) ? 1 : 0;

      GO minMaxInput[3];
      minMaxInput[0] = -minMyGID_;
      minMaxInput[1] = maxMyGID_;
      minMaxInput[2] = localDist;

      GO minMaxOutput[3];
      minMaxOutput[0] = 0;
      minMaxOutput[1] = 0;
      minMaxOutput[2] = 0;
      reduceAll<int, GO> (*comm, REDUCE_MAX, 3, minMaxInput, minMaxOutput);
      minAllGID_ = -minMaxOutput[0];
      maxAllGID_ = minMaxOutput[1];
      const GO globalDist = minMaxOutput[2];
      distributed_ = (comm_->getSize () > 1 && globalDist == 1);
    }
    else { // unsigned; use two reductions
      // This is always correct, no matter the signedness of GO.
      reduceAll<int, GO> (*comm_, REDUCE_MIN, minMyGID_, outArg (minAllGID_));
      reduceAll<int, GO> (*comm_, REDUCE_MAX, maxMyGID_, outArg (maxAllGID_));
      distributed_ = checkIsDist ();
    }

    contiguous_  = false; // "Contiguous" is conservative.

    TEUCHOS_TEST_FOR_EXCEPTION(
      minAllGID_ < indexBase_,
      std::invalid_argument,
      "Tpetra::Map constructor (noncontiguous): "
      "Minimum global ID = " << minAllGID_ << " over all process(es) is "
      "less than the given indexBase = " << indexBase_ << ".");

    // Create the Directory on demand in getRemoteIndexList().
    //setupDirectory ();
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::~Map ()
  {}


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::isOneToOne () const
  {
    // This is a collective operation, if it hasn't been called before.
    setupDirectory ();

    return directory_->isOneToOne (* (getComm ()));
  }


  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LocalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  getLocalElement(GlobalOrdinal globalIndex) const {
    if (isContiguous()) {
      if (globalIndex < getMinGlobalIndex() || globalIndex > getMaxGlobalIndex()) {
        return Teuchos::OrdinalTraits<LocalOrdinal>::invalid();
      }
      return Teuchos::as<LocalOrdinal>(globalIndex - getMinGlobalIndex());
    }
    else if (globalIndex >= firstContiguousGID_ &&
             globalIndex <= lastContiguousGID_) {
      return Teuchos::as<LocalOrdinal>(globalIndex - firstContiguousGID_);
    }
    else {
      // This returns Teuchos::OrdinalTraits<LocalOrdinal>::invalid()
      // if the given global index is not in the table.
      return glMap_->get (globalIndex);
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  GlobalOrdinal
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  getGlobalElement(LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();
    }
    if (isContiguous()) {
      return getMinGlobalIndex() + localIndex;
    }
    else {
      return lgMap_[localIndex];
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  isNodeLocalElement (LocalOrdinal localIndex) const {
    if (localIndex < getMinLocalIndex() || localIndex > getMaxLocalIndex()) {
      return false;
    } else {
      return true;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  isNodeGlobalElement (GlobalOrdinal globalIndex) const {
    return this->getLocalElement (globalIndex) !=
      Teuchos::OrdinalTraits<LocalOrdinal>::invalid ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isUniform() const {
    return uniform_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isContiguous() const {
    return contiguous_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  isCompatible (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    //
    // Tests that avoid the Boolean all-reduce below by using
    // globally consistent quantities.
    //
    if (this == &map) {
      // Pointer equality on one process always implies pointer
      // equality on all processes, since Map is immutable.
      return true;
    }
    else if (getComm ()->getSize () != map.getComm ()->getSize ()) {
      // The two communicators have different numbers of processes.
      // It's not correct to call isCompatible() in that case.  This
      // may result in the all-reduce hanging below.
      return false;
    }
    else if (getGlobalNumElements () != map.getGlobalNumElements ()) {
      // Two Maps are definitely NOT compatible if they have different
      // global numbers of indices.
      return false;
    }
    else if (isContiguous () && isUniform () &&
             map.isContiguous () && map.isUniform ()) {
      // Contiguous uniform Maps with the same number of processes in
      // their communicators, and with the same global numbers of
      // indices, are always compatible.
      return true;
    }
    else if (! isContiguous () && ! map.isContiguous () &&
             ! lgMap_.is_null () && ! map.lgMap_.is_null () &&
             lgMap_.getRawPtr () == map.lgMap_.getRawPtr ()) {
      // Noncontiguous Maps whose global index lists are nonnull and
      // have the same pointer must be the same (and therefore
      // contiguous).  (This must be the case, because Map's
      // constructor makes a deep copy of its input.  This is NOT
      // necessarily true for the Kokkos refactor version of Map!)
      //
      // Nonnull is important, since every empty list is null.  (For
      // example, consider a communicator with two processes, and two
      // Maps that share this communicator, with zero global indices
      // on the first process, and different numbers of global indices
      // on the second process.)
      return true;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(
      getGlobalNumElements () != map.getGlobalNumElements (), std::logic_error,
      "Tpetra::Map::isCompatible: There's a bug in this method.  We've already "
      "checked that this condition is true above, but it's false here.  "
      "Please report this bug to the Tpetra developers.");

    // Do both Maps have the same number of indices on each process?
    const int locallyCompat =
      (getNodeNumElements () == map.getNodeNumElements ()) ? 1 : 0;

    int globallyCompat = 0;
    reduceAll<int, int> (*comm_, REDUCE_MIN, locallyCompat, outArg (globallyCompat));
    return (globallyCompat == 1);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  locallySameAs (const Map<LocalOrdinal, GlobalOrdinal, Node>& map) const
  {
    using Teuchos::ArrayView;
    typedef GlobalOrdinal GO;
    typedef typename ArrayView<const GO>::size_type size_type;

    // If both Maps are contiguous, we can compare their GID ranges
    // easily by looking at the min and max GID on this process.
    // Otherwise, we'll compare their GID lists.  If only one Map is
    // contiguous, then we only have to call getNodeElementList() on
    // the noncontiguous Map.  (It's best to avoid calling it on a
    // contiguous Map, since it results in unnecessary storage that
    // persists for the lifetime of the Map.)

    if (getNodeNumElements () != map.getNodeNumElements ()) {
      return false;
    }
    else if (getMinGlobalIndex () != map.getMinGlobalIndex () ||
             getMaxGlobalIndex () != map.getMaxGlobalIndex ()) {
      return false;
    }
    else {
      if (isContiguous ()) {
        if (map.isContiguous ()) {
          return true; // min and max match, so the ranges match.
        }
        else { // *this is contiguous, but map is not contiguous
          TEUCHOS_TEST_FOR_EXCEPTION(
            ! this->isContiguous () || map.isContiguous (), std::logic_error,
            "Tpetra::Map::locallySameAs: BUG");
          ArrayView<const GO> rhsElts = map.getNodeElementList ();
          const GO minLhsGid = this->getMinGlobalIndex ();
          const size_type numRhsElts = rhsElts.size ();
          for (size_type k = 0; k < numRhsElts; ++k) {
            const GO curLhsGid = minLhsGid + static_cast<GO> (k);
            if (curLhsGid != rhsElts[k]) {
              return false; // stop on first mismatch
            }
          }
          return true;
        }
      }
      else if (map.isContiguous ()) { // *this is not contiguous, but map is
        TEUCHOS_TEST_FOR_EXCEPTION(
          this->isContiguous () || ! map.isContiguous (), std::logic_error,
          "Tpetra::Map::locallySameAs: BUG");
        ArrayView<const GO> lhsElts = this->getNodeElementList ();
        const GO minRhsGid = map.getMinGlobalIndex ();
        const size_type numLhsElts = lhsElts.size ();
        for (size_type k = 0; k < numLhsElts; ++k) {
          const GO curRhsGid = minRhsGid + static_cast<GO> (k);
          if (curRhsGid != lhsElts[k]) {
            return false; // stop on first mismatch
          }
        }
        return true;
      }
      else { // neither *this nor map are contiguous
        // std::equal requires that the latter range is as large as
        // the former.  We know the ranges have equal length, because
        // they have the same number of local entries.
        ArrayView<const GO> lhsElts =     getNodeElementList ();
        ArrayView<const GO> rhsElts = map.getNodeElementList ();
        return std::equal (lhsElts.begin (), lhsElts.end (), rhsElts.begin ());
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  isSameAs (const Map<LocalOrdinal,GlobalOrdinal,Node> &map) const
  {
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    //
    // Tests that avoid the Boolean all-reduce below by using
    // globally consistent quantities.
    //
    if (this == &map) {
      // Pointer equality on one process always implies pointer
      // equality on all processes, since Map is immutable.
      return true;
    }
    else if (getComm ()->getSize () != map.getComm ()->getSize ()) {
      // The two communicators have different numbers of processes.
      // It's not correct to call isSameAs() in that case.  This
      // may result in the all-reduce hanging below.
      return false;
    }
    else if (getGlobalNumElements () != map.getGlobalNumElements ()) {
      // Two Maps are definitely NOT the same if they have different
      // global numbers of indices.
      return false;
    }
    else if (getMinAllGlobalIndex () != map.getMinAllGlobalIndex () ||
             getMaxAllGlobalIndex () != map.getMaxAllGlobalIndex () ||
             getIndexBase () != map.getIndexBase ()) {
      // If the global min or max global index doesn't match, or if
      // the index base doesn't match, then the Maps aren't the same.
      return false;
    }
    else if (isDistributed () != map.isDistributed ()) {
      // One Map is distributed and the other is not, which means that
      // the Maps aren't the same.
      return false;
    }
    else if (isContiguous () && isUniform () &&
             map.isContiguous () && map.isUniform ()) {
      // Contiguous uniform Maps with the same number of processes in
      // their communicators, with the same global numbers of indices,
      // and with matching index bases and ranges, must be the same.
      return true;
    }

    // The two communicators must have the same number of processes,
    // with process ranks occurring in the same order.  This uses
    // MPI_COMM_COMPARE.  The MPI 3.1 standard (Section 6.4) says:
    // "Operations that access communicators are local and their
    // execution does not require interprocess communication."
    // However, just to be sure, I'll put this call after the above
    // tests that don't communicate.
    if (! Details::congruent (*comm_, * (map.getComm ()))) {
      return false;
    }

    // If we get this far, we need to check local properties and then
    // communicate local sameness across all processes.
    const int isSame_lcl = locallySameAs (map) ? 1 : 0;

    // Return true if and only if all processes report local sameness.
    int isSame_gbl = 0;
    reduceAll<int, int> (*comm_, REDUCE_MIN, isSame_lcl, outArg (isSame_gbl));
    return isSame_gbl == 1;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::ArrayView<const GlobalOrdinal>
  Map<LocalOrdinal,GlobalOrdinal,Node>::getNodeElementList () const
  {
    // If the local-to-global mapping doesn't exist yet, and if we
    // have local entries, then create and fill the local-to-global
    // mapping.
    if (lgMap_.is_null() && numLocalElements_ > 0) {
#ifdef HAVE_TEUCHOS_DEBUG
      // The local-to-global mapping should have been set up already
      // for a noncontiguous map.
      TEUCHOS_TEST_FOR_EXCEPTION( ! isContiguous(), std::logic_error,
        "Tpetra::Map::getNodeElementList: The local-to-global mapping (lgMap_) "
        "should have been set up already for a noncontiguous Map.  Please report"
        " this bug to the Tpetra team.");
#endif // HAVE_TEUCHOS_DEBUG

      typedef typename Teuchos::ArrayRCP<GlobalOrdinal>::size_type size_type;
      const size_type numElts = Teuchos::as<size_type> (getNodeNumElements ());
      lgMap_ = Teuchos::arcp<GlobalOrdinal> (numElts);

      if (numElts != 0) {
        // mfh 20 Feb 2013: Older compilers don't inline
        // ArrayRCP::operator* well, so we only use it in a debug build
        // (where its bounds checking is useful).
#ifdef HAVE_TEUCHOS_DEBUG
        Teuchos::ArrayRCP<GlobalOrdinal> lgMapPtr = lgMap_;
#else
        GlobalOrdinal* const lgMapPtr = lgMap_.getRawPtr ();
#endif // HAVE_TEUCHOS_DEBUG
        GlobalOrdinal gid = minMyGID_;
        for (size_type k = 0; k < numElts; ++k, ++gid) {
          lgMapPtr[k] = gid;
        }
      }
    }
    return lgMap_ ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::isDistributed() const {
    return distributed_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string Map<LocalOrdinal,GlobalOrdinal,Node>::description() const {
    using Teuchos::TypeNameTraits;
    std::ostringstream os;

    os << "Tpetra::Map: {"
       << "LocalOrdinalType: " << TypeNameTraits<LocalOrdinal>::name ()
       << ", GlobalOrdinalType: " << TypeNameTraits<GlobalOrdinal>::name ()
       << ", NodeType: " << TypeNameTraits<Node>::name ();
    if (this->getObjectLabel () != "") {
      os << ", Label: \"" << this->getObjectLabel () << "\"";
    }
    os << ", Global number of entries: " << getGlobalNumElements ()
       << ", Number of processes: " << getComm ()->getSize ()
       << ", Uniform: " << (isUniform () ? "true" : "false")
       << ", Contiguous: " << (isContiguous () ? "true" : "false")
       << ", Distributed: " << (isDistributed () ? "true" : "false")
       << "}";
    return os.str ();
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Map<LocalOrdinal,GlobalOrdinal,Node>::
  describe (Teuchos::FancyOStream &out,
            const Teuchos::EVerbosityLevel verbLevel) const
  {
    using std::endl;
    using std::setw;
    using Teuchos::ArrayView;
    using Teuchos::as;
    using Teuchos::OSTab;
    using Teuchos::toString;
    using Teuchos::TypeNameTraits;
    using Teuchos::VERB_DEFAULT;
    using Teuchos::VERB_NONE;
    using Teuchos::VERB_LOW;
    using Teuchos::VERB_MEDIUM;
    using Teuchos::VERB_HIGH;
    using Teuchos::VERB_EXTREME;
    typedef typename ArrayView<const GlobalOrdinal>::size_type size_type;

    const size_t nME = getNodeNumElements ();
    ArrayView<const GlobalOrdinal> myEntries = getNodeElementList ();
    const int myRank = comm_->getRank ();
    const int numProcs = comm_->getSize ();

    const Teuchos::EVerbosityLevel vl = (verbLevel == VERB_DEFAULT) ? VERB_LOW : verbLevel;

    size_t width = 1;
    for (size_t dec=10; dec<getGlobalNumElements(); dec *= 10) {
      ++width;
    }
    width = std::max<size_t> (width, as<size_t> (12)) + 2;

    // By convention, describe() always begins with a tab before printing.
    OSTab tab0 (out);

    if (vl == VERB_NONE) {
      // do nothing
    }
    else if (vl == VERB_LOW) {
      if (myRank == 0) {
        out << "Tpetra::Map:" << endl;
        OSTab tab1 (out);
        out << "LocalOrdinalType: " << TypeNameTraits<LocalOrdinal>::name () << endl
            << "GlobalOrdinalType: " << TypeNameTraits<GlobalOrdinal>::name () << endl
            << "NodeType: " << TypeNameTraits<Node>::name () << endl;
        if (this->getObjectLabel () != "") {
          out << "Label: \"" << this->getObjectLabel () << "\"" << endl;
        }
        out << "Global number of entries: " << getGlobalNumElements () << endl
            << "Minimum global index: " << getMinAllGlobalIndex () << endl
            << "Maximum global index: " << getMaxAllGlobalIndex () << endl
            << "Index base: " << getIndexBase () << endl
            << "Number of processes: " << getComm ()->getSize () << endl
            << "Uniform: " << (isUniform () ? "true" : "false") << endl
            << "Contiguous: " << (isContiguous () ? "true" : "false") << endl
            << "Distributed: " << (isDistributed () ? "true" : "false") << endl;
      }
    }

    if (vl >= VERB_HIGH) { // HIGH or EXTREME
      for (int p = 0; p < numProcs; ++p) {
        if (myRank == p) {
          out << "Process " << myRank << ":" << endl;
          OSTab tab1 (out);
          out << "My number of entries: " << nME << endl
              << "My minimum global index: " << getMinGlobalIndex () << endl
              << "My maximum global index: " << getMaxGlobalIndex () << endl;
          if (vl == VERB_EXTREME) {
            out << "My global indices: [";
            for (size_type k = 0; k < myEntries.size (); ++k) {
              out << myEntries[k];
              if (k + 1 < myEntries.size ()) {
                out << ", ";
              }
            }
            out << "]" << endl;
          }
          std::flush (out);
        }
        // Do a few global ops to give I/O a chance to complete
        comm_->barrier ();
        comm_->barrier ();
        comm_->barrier ();
      }
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  Map<LocalOrdinal, GlobalOrdinal, Node>::
  replaceCommWithSubset (const Teuchos::RCP<const Teuchos::Comm<int> >& newComm) const
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::OrdinalTraits;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;
    typedef global_size_t GST;
    typedef LocalOrdinal LO;
    typedef GlobalOrdinal GO;
    typedef Map<LO, GO, Node> map_type;

    // mfh 26 Mar 2013: The lazy way to do this is simply to recreate
    // the Map by calling its ordinary public constructor, using the
    // original Map's data.  This only involves O(1) all-reduces over
    // the new communicator, which in the common case only includes a
    // small number of processes.

    // Make Map compute the global number of elements.
    const GST globalNumElts = OrdinalTraits<GST>::invalid ();
    ArrayView<const GO> myElts = this->getNodeElementList ();
    RCP<Node> node = this->getNode ();

    // Create the Map to return.
    if (newComm.is_null ()) {
      return null; // my process does not participate in the new Map
    } else {
      // Map requires that the index base equal the global min GID.
      // Figuring out the global min GID requires a reduction over all
      // processes in the new communicator.  It could be that some (or
      // even all) of these processes contain zero entries.  (Recall
      // that this method, unlike removeEmptyProcesses(), may remove
      // an arbitrary subset of processes.)  We deal with this by
      // doing a min over the min GID on each process if the process
      // has more than zero entries, or the global max GID, if that
      // process has zero entries.  If no processes have any entries,
      // then the index base doesn't matter anyway.
      const GO myMinGid = (this->getNodeNumElements () == 0) ?
        this->getMaxAllGlobalIndex () : this->getMinGlobalIndex ();
      GO newIndexBase = OrdinalTraits<GO>::invalid ();
      reduceAll<int, GO> (*newComm, REDUCE_MIN, myMinGid, outArg (newIndexBase));
      return rcp (new map_type (globalNumElts, myElts, newIndexBase, newComm, node));
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >
  Map<LocalOrdinal, GlobalOrdinal, Node>::
  removeEmptyProcesses () const
  {
    using Teuchos::Comm;
    using Teuchos::null;
    using Teuchos::outArg;
    using Teuchos::RCP;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    // Create the new communicator.  split() returns a valid
    // communicator on all processes.  On processes where color == 0,
    // ignore the result.  Passing key == 0 tells MPI to order the
    // processes in the new communicator by their rank in the old
    // communicator.
    const int color = (numLocalElements_ == 0) ? 0 : 1;
    // MPI_Comm_split must be called collectively over the original
    // communicator.  We can't just call it on processes with color
    // one, even though we will ignore its result on processes with
    // color zero.
    RCP<const Comm<int> > newComm = comm_->split (color, 0);
    if (color == 0) {
      newComm = null;
    }

    // Create the Map to return.
    if (newComm.is_null ()) {
      return null; // my process does not participate in the new Map
    } else {
      // The default constructor that's useful for clone() above is
      // also useful here.
      RCP<Map> map            = rcp (new Map ());

      map->comm_              = newComm;
      map->indexBase_         = indexBase_;
      map->numGlobalElements_ = numGlobalElements_;
      map->numLocalElements_  = numLocalElements_;
      map->minMyGID_          = minMyGID_;
      map->maxMyGID_          = maxMyGID_;
      map->minAllGID_         = minAllGID_;
      map->maxAllGID_         = maxAllGID_;
      map->firstContiguousGID_= firstContiguousGID_;
      map->lastContiguousGID_ = lastContiguousGID_;

      // Uniformity and contiguity have not changed.  The directory
      // has changed, but we've taken care of that above.
      map->uniform_    = uniform_;
      map->contiguous_ = contiguous_;

      // If the original Map was NOT distributed, then the new Map
      // cannot be distributed.
      //
      // If the number of processes in the new communicator is 1, then
      // the new Map is not distributed.
      //
      // Otherwise, we have to check the new Map using an all-reduce
      // (over the new communicator).  For example, the original Map
      // may have had some processes with zero elements, and all other
      // processes with the same number of elements as in the whole
      // Map.  That Map is technically distributed, because of the
      // processes with zero elements.  Removing those processes would
      // make the new Map locally replicated.
      if (! distributed_ || newComm->getSize () == 1) {
        map->distributed_ = false;
      } else {
        const int iOwnAllGids = (numLocalElements_ == numGlobalElements_) ? 1 : 0;
        int allProcsOwnAllGids = 0;
        reduceAll<int, int> (*newComm, REDUCE_MIN, iOwnAllGids, outArg (allProcsOwnAllGids));
        map->distributed_ = (allProcsOwnAllGids == 1) ? false : true;
      }

      map->lgMap_ = lgMap_;
      map->glMap_ = glMap_;
      map->node_ = node_;

      // The Directory will be created on demand in getRemoteIndexList().
      //
      // FIXME (mfh 26 Mar 2013) It should be possible to "filter" the
      // directory more efficiently than just recreating it.  If
      // directory recreation proves a bottleneck, we can always
      // revisit this.  On the other hand, Directory creation is only
      // collective over the new, presumably much smaller
      // communicator, so it may not be worth the effort to optimize.
      map->directory_ = null;
      return map;
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  Map<LocalOrdinal,GlobalOrdinal,Node>::setupDirectory () const
  {
    using Teuchos::rcp;
    typedef Directory<LocalOrdinal,GlobalOrdinal,Node> directory_type;
    // Only create the Directory if it hasn't been created yet.
    // This is a collective operation.
    if (directory_.is_null ()) {
      directory_ = rcp (new directory_type (*this));
    }
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                    const Teuchos::ArrayView<int> & imageIDList,
                    const Teuchos::ArrayView<LocalOrdinal> & LIDList) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDList.size() > 0 && getGlobalNumElements() == 0, std::runtime_error,
      Teuchos::typeName(*this) << "::getRemoteIndexList: The Map has zero "
      "entries (globally), so you may not call this method.");
    // getRemoteIndexList must be called collectively, and Directory
    // creation is collective too, so it's OK to create the Directory
    // on demand.
    setupDirectory ();
    return directory_->getDirectoryEntries (*this, GIDList, imageIDList, LIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  LookupStatus Map<LocalOrdinal,GlobalOrdinal,Node>::getRemoteIndexList(
                    const Teuchos::ArrayView<const GlobalOrdinal> & GIDList,
                    const Teuchos::ArrayView<int> & imageIDList) const {
    TEUCHOS_TEST_FOR_EXCEPTION(
      GIDList.size() > 0 && getGlobalNumElements() == 0, std::runtime_error,
      Teuchos::typeName(*this) << "::getRemoteIndexList: The Map has zero "
      "entries (globally), so you may not call this method.");
    // getRemoteIndexList must be called collectively, and Directory
    // creation is collective too, so it's OK to create the Directory
    // on demand.
    setupDirectory ();
    return directory_->getDirectoryEntries (*this, GIDList, imageIDList);
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<const Teuchos::Comm<int> >
  Map<LocalOrdinal,GlobalOrdinal,Node>::getComm() const {
    return comm_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Node>
  Map<LocalOrdinal,GlobalOrdinal,Node>::getNode() const {
    return node_;
  }

  template <class LocalOrdinal,class GlobalOrdinal, class Node>
  bool Map<LocalOrdinal,GlobalOrdinal,Node>::checkIsDist() const {
    using Teuchos::as;
    using Teuchos::outArg;
    using Teuchos::REDUCE_MIN;
    using Teuchos::reduceAll;

    bool global = false;
    if (comm_->getSize () > 1) {
      // The communicator has more than one process, but that doesn't
      // necessarily mean the Map is distributed.
      int localRep = 0;
      if (numGlobalElements_ == as<global_size_t> (numLocalElements_)) {
        // The number of local elements on this process equals the
        // number of global elements.
        //
        // NOTE (mfh 22 Nov 2011) Does this still work if there were
        // duplicates in the global ID list on input (the third Map
        // constructor), so that the number of local elements (which
        // are not duplicated) on this process could be less than the
        // number of global elements, even if this process owns all
        // the elements?
        localRep = 1;
      }
      int allLocalRep;
      reduceAll<int, int> (*comm_, REDUCE_MIN, localRep, outArg (allLocalRep));
      if (allLocalRep != 1) {
        // At least one process does not own all the elements.
        // This makes the Map a distributed Map.
        global = true;
      }
    }
    // If the communicator has only one process, then the Map is not
    // distributed.
    return global;
  }

} // Tpetra namespace

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> >
Tpetra::createLocalMap(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return createLocalMapWithNode<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType>(numElements, comm, KokkosClassic::DefaultNode::getDefaultNode());
}

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> >
Tpetra::createUniformContigMap(global_size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return createUniformContigMapWithNode<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType>(numElements, comm, KokkosClassic::DefaultNode::getDefaultNode());
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createUniformContigMapWithNode (global_size_t numElements,
                                        const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                        const Teuchos::RCP<Node>& node)
{
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  const GlobalOrdinal indexBase = Teuchos::OrdinalTraits<GlobalOrdinal>::zero ();

  return rcp (new map_type (numElements, indexBase, comm, GloballyDistributed, node));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createLocalMapWithNode(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  using Tpetra::global_size_t;
  using Teuchos::as;
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  const GlobalOrdinal indexBase = Teuchos::OrdinalTraits<GlobalOrdinal>::zero ();
  const global_size_t globalNumElts = as<global_size_t> (numElements);

  return rcp (new map_type (globalNumElts, indexBase, comm, LocallyReplicated, node));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createContigMapWithNode(Tpetra::global_size_t numElements, size_t localNumElements,
                                const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  const GlobalOrdinal indexBase = Teuchos::OrdinalTraits<GlobalOrdinal>::zero ();

  return rcp (new map_type (numElements, localNumElements, indexBase, comm, node));
}

template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> >
Tpetra::createContigMap(Tpetra::global_size_t numElements, size_t localNumElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm) {
  return Tpetra::createContigMapWithNode<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType>(numElements, localNumElements, comm, KokkosClassic::DefaultNode::getDefaultNode() );
}


template <class LocalOrdinal, class GlobalOrdinal>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType> >
Tpetra::createNonContigMap(const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
                           const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  return Tpetra::createNonContigMapWithNode<LocalOrdinal,GlobalOrdinal,KokkosClassic::DefaultNode::DefaultNodeType>(elementList, comm, KokkosClassic::DefaultNode::getDefaultNode() );
}


template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createNonContigMapWithNode(const Teuchos::ArrayView<const GlobalOrdinal> &elementList,
                                   const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                                   const Teuchos::RCP<Node> &node)
{
  using Teuchos::rcp;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> map_type;
  typedef Tpetra::global_size_t GST;
  return rcp (new map_type (Teuchos::OrdinalTraits<GST>::invalid (),
                            elementList,
                            Teuchos::OrdinalTraits<GST>::zero (),
                            comm,
                            node));
}

template <class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP< const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createWeightedContigMapWithNode(int myWeight, Tpetra::global_size_t numElements,
                                        const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< Node > &node) {
  Teuchos::RCP< Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > map;
  int sumOfWeights, elemsLeft, localNumElements;
  const int numImages = comm->getSize(),
            myImageID = comm->getRank();
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,myWeight,Teuchos::outArg(sumOfWeights));
  const double myShare = ((double)myWeight) / ((double)sumOfWeights);
  localNumElements = (int)std::floor( myShare * ((double)numElements) );
  // std::cout << "numElements: " << numElements << "  myWeight: " << myWeight << "  sumOfWeights: " << sumOfWeights << "  myShare: " << myShare << std::endl;
  Teuchos::reduceAll<int>(*comm,Teuchos::REDUCE_SUM,localNumElements,Teuchos::outArg(elemsLeft));
  elemsLeft = numElements - elemsLeft;
  // std::cout << "(before) localNumElements: " << localNumElements << "  elemsLeft: " << elemsLeft << std::endl;
  // i think this is true. just test it for now.
  TEUCHOS_TEST_FOR_EXCEPT(elemsLeft < -numImages || numImages < elemsLeft);
  if (elemsLeft < 0) {
    // last elemsLeft nodes lose an element
    if (myImageID >= numImages-elemsLeft) --localNumElements;
  }
  else if (elemsLeft > 0) {
    // first elemsLeft nodes gain an element
    if (myImageID < elemsLeft) ++localNumElements;
  }
  // std::cout << "(after) localNumElements: " << localNumElements << std::endl;
  return createContigMapWithNode<LocalOrdinal,GlobalOrdinal,Node>(numElements,localNumElements,comm,node);

}


template<class LO, class GO, class NT>
Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >
Tpetra::createOneToOne (const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& M)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::rcp;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef global_size_t GST;
  const GST GINV = Teuchos::OrdinalTraits<GST>::invalid ();
  const int myRank = M->getComm ()->getRank ();

  // Bypasses for special cases where either M is known to be
  // one-to-one, or the one-to-one version of M is easy to compute.
  // This is why we take M as an RCP, not as a const reference -- so
  // that we can return M itself if it is 1-to-1.
  if (! M->isDistributed ()) {
    // For a locally replicated Map, we assume that users want to push
    // all the GIDs to Process 0.

    // mfh 05 Nov 2013: getGlobalNumElements() does indeed return what
    // you think it should return, in this special case of a locally
    // replicated contiguous Map.
    const GST numGlobalEntries = M->getGlobalNumElements ();
    if (M->isContiguous ()) {
      const size_t numLocalEntries =
        (myRank == 0) ? as<size_t> (numGlobalEntries) : static_cast<size_t> (0);
      return rcp (new map_type (numGlobalEntries, numLocalEntries,
                                M->getIndexBase (), M->getComm (),
                                M->getNode ()));
    }
    else {
      ArrayView<const GO> myGids =
        (myRank == 0) ? M->getNodeElementList () : Teuchos::null;
      return rcp (new map_type (GINV, myGids (), M->getIndexBase (),
                                M->getComm (), M->getNode ()));

    }
  }
  else if (M->isContiguous ()) {
    // Contiguous, distributed Maps are one-to-one by construction.
    // (Locally replicated Maps can be contiguous.)
    return M;
  }
  else {
    Tpetra::Directory<LO, GO, NT> directory (*M);
    const size_t numMyElems = M->getNodeNumElements ();
    ArrayView<const GO> myElems = M->getNodeElementList ();
    Array<int> owner_procs_vec (numMyElems);

    directory.getDirectoryEntries (*M, myElems, owner_procs_vec ());

    Array<GO> myOwned_vec (numMyElems);
    size_t numMyOwnedElems = 0;
    for (size_t i = 0; i < numMyElems; ++i) {
      const GO GID = myElems[i];
      const int owner = owner_procs_vec[i];

      if (myRank == owner) {
        myOwned_vec[numMyOwnedElems++] = GID;
      }
    }
    myOwned_vec.resize (numMyOwnedElems);

    return rcp (new map_type (GINV, myOwned_vec (), M->getIndexBase (),
                              M->getComm (), M->getNode ()));
  }
}

template<class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
Tpetra::createOneToOne (const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > &M,
                        const Tpetra::Details::TieBreak<LocalOrdinal,GlobalOrdinal> & tie_break)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::rcp;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO,GO,Node> map_type;
  int myID = M->getComm()->getRank();

  // FIXME (mfh 20 Feb 2013) We should have a bypass for contiguous
  // Maps (which are 1-to-1 by construction).

  //Based off Epetra's one to one.

  Tpetra::Directory<LO, GO, Node> directory (*M, tie_break);
  size_t numMyElems = M->getNodeNumElements ();
  ArrayView<const GO> myElems = M->getNodeElementList ();
  Array<int> owner_procs_vec (numMyElems);

  directory.getDirectoryEntries (*M, myElems, owner_procs_vec ());

  Array<GO> myOwned_vec (numMyElems);
  size_t numMyOwnedElems = 0;
  for (size_t i = 0; i < numMyElems; ++i) {
    GO GID = myElems[i];
    int owner = owner_procs_vec[i];

    if (myID == owner) {
      myOwned_vec[numMyOwnedElems++] = GID;
    }
  }
  myOwned_vec.resize (numMyOwnedElems);

  const global_size_t GINV =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();
  return rcp (new map_type (GINV, myOwned_vec (), M->getIndexBase (),
                            M->getComm (), M->getNode ()));
}

#if defined(TPETRA_USE_KOKKOS_REFACTOR_MAP)
// Include KokkosRefactor partial specialization if enabled
#  if defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#    include "Tpetra_KokkosRefactor_Map_def.hpp"
#  endif // defined(TPETRA_HAVE_KOKKOS_REFACTOR)
#endif // defined(TPETRA_USE_KOKKOS_REFACTOR_MAP)

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

//! Explicit instantiation macro supporting the Map class. Instantiates the class and the non-member constructors.
#define TPETRA_MAP_INSTANT(LO,GO,NODE) \
  \
  template class Map< LO , GO , NODE >; \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createLocalMapWithNode<LO,GO,NODE>(size_t numElements, const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createContigMapWithNode<LO,GO,NODE>(global_size_t numElements, size_t localNumElements, \
                                      const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createNonContigMapWithNode(const Teuchos::ArrayView<const GO> &elementList, \
                             const RCP<const Teuchos::Comm<int> > &comm,      \
                             const RCP<NODE> &node);                          \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createUniformContigMapWithNode<LO,GO,NODE>(global_size_t numElements, \
                                             const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP< const Map<LO,GO,NODE> > \
  createWeightedContigMapWithNode<LO,GO,NODE>(int thisNodeWeight, global_size_t numElements, \
                                              const Teuchos::RCP< const Teuchos::Comm< int > > &comm, const Teuchos::RCP< NODE > &node); \
  \
  template Teuchos::RCP<const Map<LO,GO,NODE> > \
  createOneToOne (const Teuchos::RCP<const Map<LO,GO,NODE> > &M); \
  \
  template Teuchos::RCP<const Map<LO,GO,NODE> > \
  createOneToOne (const Teuchos::RCP<const Map<LO,GO,NODE> > &M, \
                  const Tpetra::Details::TieBreak<LO,GO> & tie_break);


//! Explicit instantiation macro supporting the Map class, on the default node for specified ordinals.
#define TPETRA_MAP_INSTANT_DEFAULTNODE(LO,GO) \
  template Teuchos::RCP< const Map<LO,GO> > \
  createLocalMap<LO,GO>( size_t, const Teuchos::RCP< const Teuchos::Comm< int > > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> > \
  createContigMap<LO,GO>( global_size_t, size_t, \
                          const Teuchos::RCP< const Teuchos::Comm< int > > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> >  \
  createNonContigMap(const Teuchos::ArrayView<const GO> &,          \
                     const RCP<const Teuchos::Comm<int> > &); \
  \
  template Teuchos::RCP< const Map<LO,GO> >  \
  createUniformContigMap<LO,GO>( global_size_t, \
                                 const Teuchos::RCP< const Teuchos::Comm< int > > &); \

#endif // TPETRA_MAP_DEF_HPP
