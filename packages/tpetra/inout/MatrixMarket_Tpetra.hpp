//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

#ifndef __MatrixMarket_Tpetra_hpp
#define __MatrixMarket_Tpetra_hpp

/// \file MatrixMarket_Tpetra.hpp
/// \brief Matrix Market file reader for Tpetra::CrsMatrix.
///

#include "Tpetra_CrsMatrix.hpp"
#include "MatrixMarket_raw.hpp"
#include "MatrixMarket_Banner.hpp"
#include "MatrixMarket_CoordDataReader.hpp"
#include "MatrixMarket_util.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <stdexcept>

namespace Tpetra {
  namespace MatrixMarket {

    /// \class Reader
    /// \brief Matrix Market file reader for Tpetra::CrsMatrix.
    ///
    /// The readSparse() and readSparseFile() class methods read in a
    /// Matrix Market "coordinate" format sparse matrix input stream
    /// resp. file (valid on MPI Rank 0), and return a
    /// Tpetra::CrsMatrix sparse matrix (SparseMatrixType) distributed
    /// in block row fashion over all the MPI ranks represented in the
    /// given communicator.
    template<class SparseMatrixType>
    class Reader {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;
      /// \typedef scalar_type
      /// \brief Type of the entries of the sparse matrix.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// \typedef local_ordinal_type
      /// \brief Only used to define map_type.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// \typedef global_ordinal_type
      /// \brief Type of indices as read from the Matrix Market file.
      ///
      /// Indices of the sparse matrix are read in as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      typedef typename SparseMatrixType::node_type node_type;

      typedef Teuchos::RCP<node_type> node_ptr;
      typedef Teuchos::Comm<int> comm_type;
      typedef Teuchos::RCP<const comm_type> comm_ptr;
      typedef Map<local_ordinal_type, global_ordinal_type, node_type> map_type;
      typedef Teuchos::RCP<const map_type> map_ptr;

    private:

      /// \typedef size_type
      /// \brief Handy typedef for entries of arrays such as rowPtr.
      typedef typename ArrayRCP<global_ordinal_type>::size_type size_type;

      /// \brief Compute initial range map.
      ///
      /// Range maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      ///
      /// \param pComm [in] Global communicator.
      /// \param pNode [in] Kokkos Node object.
      /// \param numRows [in] Global number of rows in the matrix.
      ///
      /// \return Range map to be used for constructing a CrsMatrix.
      ///
      static Teuchos::RCP<const map_type>
      makeRangeMap (const Teuchos::RCP<const comm_type>& pComm,
		    const Teuchos::RCP<node_type>& pNode,
		    const global_ordinal_type numRows)
      {
	using Teuchos::rcp;

	// A conventional, uniformly partitioned, contiguous map.
	return rcp (new map_type (static_cast<global_size_t> (numRows), 
				  static_cast<global_ordinal_type> (0),
				  pComm, GloballyDistributed, pNode));
      }

      /// \brief Compute initial row map, or verify an existing one.
      ///
      /// The typical case when reading a sparse matrix from a file is
      /// for the reader itself to create a new row map, in particular
      /// a standard uniform contiguous one-to-one row map.  However,
      /// we also give the option to use an existing row map, if you
      /// are already using a particular distribution for (say) vector
      /// data and don't want to stop using it.  In the latter case
      /// (pRowMap is not null), we validate the communicator and node
      /// of the existing row map that you pass in.  In either case,
      /// you need to know the (global) number of rows in the matrix.
      ///
      /// \param pRowMap [in] If non-null, test pRowMap for validity,
      ///   and return it if valid.  Otherwise, if pRowMap is null,
      ///   initialize and return a (uniform contiguous one-to-one)
      ///   row map.  "Validity" here means that the map's
      ///   communicator and node are the same objects (pointerwise)
      ///   as the corresponding arguments.  (Note that the global
      ///   number of elements may not be the same as the number of
      ///   rows; a row map is not required to be one-to-one.)  The
      ///   typical case is to pass in null here, which is why we call
      ///   this routine "makeRowMap".
      ///
      /// \param pComm [in] Global communicator.
      ///
      /// \param pNode [in] Kokkos Node object.
      ///
      /// \param numRows [in] Global number of rows in the matrix.  If
      ///   pRowMap is nonnull, used only for error checking.
      ///
      /// \return If pRowMap is null, a new row map, otherwise pRowMap.
      ///
      static Teuchos::RCP<const map_type>
      makeRowMap (const Teuchos::RCP<const map_type>& pRowMap, 
		  const Teuchos::RCP<const comm_type>& pComm,
		  const Teuchos::RCP<node_type>& pNode,
		  const global_ordinal_type numRows)
      {
	using Teuchos::rcp;

	// If the caller didn't provide a map, return a conventional,
	// uniformly partitioned, contiguous map.
	if (pRowMap.is_null())
	  return rcp (new map_type (static_cast<global_size_t> (numRows), 
				    static_cast<global_ordinal_type> (0),
				    pComm, GloballyDistributed, pNode));
	else 
	  {
	    TEST_FOR_EXCEPTION(! pRowMap->isDistributed() && pComm->getSize() > 1, 
			       std::invalid_argument,
			       "The specified row map is not distributed, but the "
			       "given communicator includes > 1 ranks.");
	    TEST_FOR_EXCEPTION(pRowMap->getComm() != pComm, 
			       std::invalid_argument,
			       "The specified row map's communicator (pRowMap->"
			       "getComm()) is different than the given separately "
			       "supplied communicator pComm.");
	    TEST_FOR_EXCEPTION(pRowMap->getNode() != pNode,
			       std::invalid_argument,
			       "The specified row map's node (pRowMap->getNode()) "
			       "is different than the given separately supplied "
			       "node pNode.");
	    return pRowMap;
	  }
      }

      /// \brief Compute domain map.
      ///
      /// Domain maps must always be one-to-one.  We will use this map
      /// when we call fillComplete() on the CrsMatrix that the reader
      /// constructs.
      /// 
      /// \param pRangeMap [in] Valid range map of the matrix, 
      ///   as returned by \c makeRangeMap().
      /// \param numRows [in] Global number of rows in the matrix.
      /// \param numCols [in] Global number of columns in the matrix.
      ///
      /// \return The domain map.  If numRows == numCols, this is
      ///   identical to the range map, otherwise we make a new map
      ///   for the domain.
      static map_ptr
      makeDomainMap (const map_ptr& pRangeMap,
		     const global_ordinal_type numRows,
		     const global_ordinal_type numCols)
      {
	typedef local_ordinal_type LO;
	typedef global_ordinal_type GO;
	typedef node_type Node;

	if (numRows == numCols) 
	  return pRangeMap;
	else
	  {
	    comm_ptr pComm = pRangeMap->getComm();
	    node_ptr pNode = pRangeMap->getNode();
	    return createUniformContigMapWithNode<LO,GO,Node> (numCols, pComm, pNode);
	  }
      }

      /// \brief Distribute the sparse matrix entries.
      ///
      /// This is one of those routines that just has to be messy.  We try
      /// to comment everything so you can see what's going on.  All of the
      /// matrix data starts out on Rank 0, stored in CSR format (rowPtr,
      /// colInd, values).  numEntriesPerRow applies to this data; it can be
      /// computed from rowPtr, and has one less entry than rowPtr.  Rank 0
      /// then engages in a dialog with Ranks 1 .. numProcs-1 to determine
      /// which part of the matrix data is theirs, and to send them their
      /// data.
      ///
      /// \param myNumEntriesPerRow [out] For my row indices, the number of
      ///   entries per row.  This array has
      ///   pRowMap->getNodeElementList().size() entries, and is indexed in
      ///   the same way, so that myNumEntriesPerRow[k] is the number of
      ///   entries in row pRowMap->getNodeElementList()[k].
      ///
      /// \param myRowPtr [out] The row pointer array for the rows
      ///   that belong to me.  This array has one more entry than
      ///   myNumEntriesPerRow, and is the prefix sum (with
      ///   myRowPtr[0] = 0) of myNumEntriesPerRow.
      ///
      /// \param myColInd [out] My rows' column indices.  If myRows =
      ///   pRowMap->getNodeElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then 
      ///   myColInd[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param myValues [out] My rows' stored matrix values.  If myRows =
      ///   pRowMap->getNodeElementList(), start = sum(myNumEntriesPerRow[0
      ///   .. k-1]), and end = start + myNumEntriesPerRow[k], then 
      ///   myValues[start .. end-1] are the column indices for myRows[k].
      ///
      /// \param pRowMap [in] Map describing the distribution of rows among
      ///   processors.
      ///
      /// \param numEntriesPerRow [in/out] For all row indices, the number
      ///   of entries per row.  You can construct this from the usual CSR
      ///   matrix "rowPtr" array by differences: numEntriesPerRow[k] =
      ///   rowPtr[k+1] - rowPtr[k].  As a result, numEntriesPerRow has one
      ///   less entry than rowPtr.  On output, the reference is invalidated
      ///   to save space.
      ///
      /// \param rowPtr [in/out] On input: the usual CSR matrix row pointer
      ///   array for the whole matrix.  On output, the reference is
      ///   invalidated to save space.
      ///
      /// \param colInd [in/out] On input: all rows' column indices.  If
      ///   start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then colInd[start .. end-1] are the column
      ///   indices for row k.  On output, the reference is invalidated to
      ///   save space.
      ///
      /// \param values [in/out] On input: all rows' stored matrix values.
      ///   If start = sum(numEntriesPerRow[0 .. k-1]), and end = start +
      ///   numEntriesPerRow[k], then values[start .. end-1] are the values
      ///   for row k.  On output, the reference is invalidated to save
      ///   space.
      ///
      /// \param debug [in] If true, print copious debugging output to
      ///   stderr on Rank 0.  This option is unlikely to be useful to
      ///   anyone but a Tpetra developer debugging this code.
      /// 
      static void
      distribute (ArrayRCP<size_t>& myNumEntriesPerRow,
		  ArrayRCP<size_type>& myRowPtr,
		  ArrayRCP<global_ordinal_type>& myColInd,
		  ArrayRCP<scalar_type>& myValues,
		  const Teuchos::RCP<const map_type>& pRowMap,
		  ArrayRCP<size_t>& numEntriesPerRow,
		  ArrayRCP<size_type>& rowPtr,
		  ArrayRCP<global_ordinal_type>& colInd,
		  ArrayRCP<scalar_type>& values,
		  const bool debug=false)
      {
	 using Teuchos::arcp;
	 using Teuchos::ArrayRCP;
	 using Teuchos::Array;
	 using Teuchos::ArrayView;
	 using Teuchos::CommRequest;
	 using Teuchos::Comm;
	 using Teuchos::RCP;
	 using Teuchos::rcpFromRef;
	 using Teuchos::receive;
	 using Teuchos::send;
	 using std::cerr;
	 using std::endl;

	 const bool extraDebug = false;
	 comm_ptr pComm = pRowMap->getComm();
	 const int numProcs = Teuchos::size (*pComm);
	 const int myRank = Teuchos::rank (*pComm);
	 const int rootRank = 0;

	 // List of the global indices of my rows.  They may or may
	 // not be contiguous, and the row map need not be one-to-one.
	 ArrayView<const global_ordinal_type> myRows = 
	   pRowMap->getNodeElementList();
	 const size_type myNumRows = myRows.size();
	 TEST_FOR_EXCEPTION(static_cast<size_t>(myNumRows) != pRowMap->getNodeNumElements(),
			    std::logic_error,
			    "pRowMap->getNodeElementList().size() = " 
			    << myNumRows
			    << " != pRowMap->getNodeNumElements() = "
			    << pRowMap->getNodeNumElements() << ".  "
			    "Please report this bug to the Tpetra developers.");
	 TEST_FOR_EXCEPTION(numEntriesPerRow.size() < myNumRows,
			    std::logic_error,
			    "numEntriesPerRow.size() = "
			    << numEntriesPerRow.size()
			    << " != pRowMap->getNodeElementList().size() = " 
			    << myNumRows << ".  Please report this bug to the "
			    "Tpetra developers.");

	 // Space for my proc's number of entries per row.  
	 // Will be filled in below.
	 myNumEntriesPerRow = arcp<size_t> (myNumRows);

	 // Teuchos::receive() returns an int; here is space for it.
	 int recvResult = 0;

	 if (myRank != rootRank)
	   {
	     // Tell the root how many rows we have.  If we're sending none,
	     // then we don't have anything else to send, nor does the root
	     // have to receive anything else.
	     send (*pComm, myNumRows, rootRank);
	     if (myNumRows != 0)
	       {
		 // Now send my rows' global indices.  Hopefully the
		 // cast to int doesn't overflow.  This is unlikely,
		 // since it should fit in a local_ordinal_type, even
		 // though it is a global_ordinal_type.
		 send (*pComm, static_cast<int> (myNumRows), 
		       myRows.getRawPtr(), rootRank);

		 // I (this proc) don't care if my global row indices are
		 // contiguous, though the root proc does (since otherwise it
		 // needs to pack noncontiguous data into contiguous storage
		 // before sending).  That's why we don't check for
		 // contiguousness here.

		 // Ask the root processor for my part of the array of the
		 // number of entries per row.
		 recvResult = receive (*pComm, rootRank, 
				       static_cast<int> (myNumRows),
				       myNumEntriesPerRow.getRawPtr());

		 // Use the resulting array to figure out how many column
		 // indices and values for which I should ask from the root
		 // processor.  
		 const size_t myNumEntries = 
		   std::accumulate (myNumEntriesPerRow.begin(), 
				    myNumEntriesPerRow.end(), 
				    static_cast<size_t> (0));

		 // Make space for my entries of the sparse matrix.  Note that
		 // they don't have to be sorted by row index.  Iterating through
		 // all my rows requires computing a running sum over
		 // myNumEntriesPerRow.
		 myColInd = arcp<global_ordinal_type> (myNumEntries);
		 myValues = arcp<scalar_type> (myNumEntries);
		 if (myNumEntries > 0)
		   { // Ask for that many column indices and values, if there are any.
		     recvResult = receive (*pComm, rootRank, 
					   static_cast<int> (myNumEntries), 
					   myColInd.getRawPtr());
		     recvResult = receive (*pComm, rootRank, 
					   static_cast<int> (myNumEntries), 
					   myValues.getRawPtr());
		   }
	       } // If I own at least one row
	   } // If I am not the root processor
	 else // I _am_ the root processor
	   {
	     if (debug)
	       cerr << "-- Proc 0: Copying my data from global arrays" << endl;

	     // Proc 0 (the root processor) still needs to (allocate,
	     // if not done already) and fill its part of the matrix
	     // (my*).
	     for (size_type k = 0; k < myNumRows; ++k)
	       {
		 const global_ordinal_type myCurRow = myRows[k];
		 const size_t numEntriesInThisRow = numEntriesPerRow[myCurRow];
		 //myNumEntriesPerRow[k] = numEntriesPerRow[myCurRow];
		 myNumEntriesPerRow[k] = numEntriesInThisRow;
	       }
	     if (extraDebug && debug)
	       {
		 cerr << "Proc " << Teuchos::rank (*(pRowMap->getComm())) 
		      << ": myNumEntriesPerRow[0.." << (myNumRows-1) << "] = [";
		 for (size_type k = 0; k < myNumRows; ++k)
		   {
		     cerr << myNumEntriesPerRow[k];
		     if (k < myNumRows-1)
		       cerr << " ";
		   }
		 cerr << "]" << endl;
	       }
	     // The total number of matrix entries that my proc owns.
	     const size_t myNumEntries = 
	       std::accumulate (myNumEntriesPerRow.begin(), 
				myNumEntriesPerRow.end(), 
				static_cast<size_t> (0));
	     if (debug)
	       cerr << "-- Proc 0: I own " << myNumRows << " rows and " 
		    << myNumEntries << " entries" << endl;

	     myColInd = arcp<global_ordinal_type> (myNumEntries);
	     myValues = arcp<scalar_type> (myNumEntries);

	     // Copy Proc 0's part of the matrix into the my* arrays.
	     // It's important that myCurPos be updated _before_ k,
	     // otherwise myCurPos will get the wrong number of entries
	     // per row (it should be for the row in the just-completed
	     // iteration, not for the next iteration's row).
	     size_t myCurPos = 0;
	     for (size_type k = 0; k < myNumRows; 
		  myCurPos += myNumEntriesPerRow[k], ++k)
	       {
		 const size_t curNumEntries = myNumEntriesPerRow[k];
		 const global_ordinal_type myRow = myRows[k];
		 const size_t curPos = rowPtr[myRow];
		 if (extraDebug && debug)
		   {
		     cerr << "k = " << k << ", myRow = " << myRow << ": colInd(" 
			  << curPos << "," << curNumEntries << "), myColInd(" 
			  << myCurPos << "," << curNumEntries << ")" << endl;
		   }
		 // Only copy if there are entries to copy, in order
		 // not to construct empty ranges for the ArrayRCP
		 // views.
		 if (curNumEntries > 0)
		   {
		     ArrayView<global_ordinal_type> colIndView = 
		       colInd(curPos, curNumEntries);
		     ArrayView<global_ordinal_type> myColIndView = 
		       myColInd(myCurPos, curNumEntries);
		     std::copy (colIndView.begin(), colIndView.end(), 
				myColIndView.begin());

		     ArrayView<scalar_type> valuesView = 
		       values(curPos, curNumEntries);
		     ArrayView<scalar_type> myValuesView = 
		       myValues(myCurPos, curNumEntries);
		     std::copy (valuesView.begin(), valuesView.end(), 
				myValuesView.begin());
		   }
	       }

	     // Proc 0 processes each other proc p in turn.
	     for (int p = 1; p < numProcs; ++p)
	       {
		 if (debug)
		   cerr << "-- Proc 0: Processing proc " << p << endl;

		 size_type theirNumRows = 0;
		 // Ask Proc p how many rows it has.  If it doesn't have any,
		 // we can move on to the next proc.  This has to be a
		 // standard receive so that we can avoid the degenerate case
		 // of sending zero data.
		 recvResult = receive (*pComm, p, &theirNumRows);
		 if (debug)
		   cerr << "-- Proc 0: Proc " << p << " owns " 
			<< theirNumRows << " rows" << endl;
		 if (theirNumRows != 0)
		   {
		     // Ask Proc p which rows it owns.  The resulting
		     // global row indices are not guaranteed to be
		     // contiguous or sorted.  Global row indices are
		     // themselves indices into the numEntriesPerRow
		     // array.

		     ArrayRCP<global_ordinal_type> theirRows = 
		       arcp<global_ordinal_type> (theirNumRows);
		     recvResult = receive (*pComm, p, 
					   static_cast<int> (theirNumRows), 
					   theirRows.getRawPtr());
		     // Extra test to make sure that the rows we
		     // received are all sensible.  This is a good
		     // idea since we are going to use the global row
		     // indices we've received to index into the
		     // numEntriesPerRow array.  Better to catch any
		     // bugs here and print a sensible error message,
		     // rather than segfault and print a cryptic error
		     // message.
		     {
		       const global_size_t numRows = 
			 pRowMap->getGlobalNumElements();
		       bool theirRowsValid = true;
		       for (size_type k = 0; k < theirNumRows; ++k)
			 {
			   // global_ordinal_t is generally a signed type.
			   if (theirRows[k] < 0)
			     theirRowsValid = false;
			   // Signed to unsigned cast should not
			   // overflow.  We cast in order to avoid
			   // compiler warnings about comparing signed
			   // and unsigned types, in case
			   // global_size_t is unsigned.
			   else if (static_cast<global_size_t> (theirRows[k]) >= numRows)
			     theirRowsValid = false;
			 }
		       if (! theirRowsValid)
			 {
			   std::ostringstream os;
			   std::copy (theirRows.begin(), theirRows.end(), 
				      std::ostream_iterator<global_ordinal_type>(os, " "));
			   TEST_FOR_EXCEPTION(! theirRowsValid, 
					      std::logic_error,
					      "Proc " << p << " has at least "
					      "one invalid row index.  Here are"
					      " all of them: [ " << os.str() 
					      << "]");
			 }
		     }

		     // Perhaps we could save a little work if we
		     // check whether Proc p's row indices are
		     // contiguous.  That would make lookups in the
		     // global data arrays faster.  For now, we just
		     // implement the general case and don't
		     // prematurely optimize.  (Remember that you're
		     // making Proc 0 read the whole file, so you've
		     // already lost scalability.)

		     // Compute the number of entries in each of Proc
		     // p's rows.  (Proc p will compute its row
		     // pointer array on its own, after it gets the
		     // data from Proc 0.)
		     ArrayRCP<size_t> theirNumEntriesPerRow;
		     theirNumEntriesPerRow = arcp<size_t> (theirNumRows);
		     for (size_type k = 0; k < theirNumRows; ++k)
		       theirNumEntriesPerRow[k] = numEntriesPerRow[theirRows[k]];

		     // Tell Proc p the number of entries in each of
		     // its rows.  Hopefully the cast to int doesn't
		     // overflow.  This is unlikely, since it should
		     // fit in a local_ordinal_type, even though it is
		     // a global_ordinal_type.
		     send (*pComm, static_cast<int> (theirNumRows), 
			   theirNumEntriesPerRow.getRawPtr(), p);

		     // Figure out how many entries Proc p owns.
		     const size_t theirNumEntries = 
		       std::accumulate (theirNumEntriesPerRow.begin(),
					theirNumEntriesPerRow.end(), 
					static_cast<size_t> (0));

		     if (debug)
		       cerr << "-- Proc 0: Proc " << p << " owns " 
			    << theirNumEntries << " entries" << endl;

		     // If there are no entries to send, then we're
		     // done with Proc p.
		     if (theirNumEntries == 0)
		       continue; 

		     // Construct (views of) proc p's column indices
		     // and values.  Later, we might like to optimize
		     // for the (common) contiguous case, for which we
		     // don't need to copy data into separate "their*"
		     // arrays (we can just use contiguous views of
		     // the global arrays).
		     ArrayRCP<global_ordinal_type> theirColInd = 
		       arcp<global_ordinal_type> (theirNumEntries);
		     ArrayRCP<scalar_type> theirValues = 
		       arcp<scalar_type> (theirNumEntries);
		     // Copy Proc p's part of the matrix into the
		     // their* arrays.  It's important that
		     // theirCurPos be updated _before_ k, otherwise
		     // theirCurPos will get the wrong number of
		     // entries per row (it should be for the row in
		     // the just-completed iteration, not for the next
		     // iteration's row).
		     size_t theirCurPos = 0;
		     for (size_type k = 0; k < theirNumRows; 
			  theirCurPos += theirNumEntriesPerRow[k], k++)
		       {
			 const size_t curNumEntries = theirNumEntriesPerRow[k];
			 const global_ordinal_type theirRow = theirRows[k];
			 const size_t curPos = rowPtr[theirRow];

			 // Only copy if there are entries to copy, in
			 // order not to construct empty ranges for
			 // the ArrayRCP views.
			 if (curNumEntries > 0)
			   {
			     ArrayView<global_ordinal_type> colIndView = 
			       colInd(curPos, curNumEntries);
			     ArrayView<global_ordinal_type> theirColIndView = 
			       theirColInd(theirCurPos, curNumEntries);
			     std::copy (colIndView.begin(), colIndView.end(), 
					theirColIndView.begin());
			 
			     ArrayView<scalar_type> valuesView = 
			       values(curPos, curNumEntries);
			     ArrayView<scalar_type> theirValuesView = 
			       theirValues(theirCurPos, curNumEntries);
			     std::copy (valuesView.begin(), valuesView.end(), 
					theirValuesView.begin());
			   }
		       }
		     // Send Proc p its column indices and values.
		     // Hopefully the cast to int doesn't overflow.
		     // This is unlikely, since it should fit in a
		     // local_ordinal_type, even though it is a
		     // global_ordinal_type.
		     send (*pComm, static_cast<int> (theirNumEntries), 
			   theirColInd.getRawPtr(), p);
		     send (*pComm, static_cast<int> (theirNumEntries), 
			   theirValues.getRawPtr(), p);

		     if (debug)
		       cerr << "-- Proc 0: Finished with proc " << p << endl;
		   } // If proc p owns at least one row
	       } // For each proc p not the root proc 0
	   } // If I'm (not) the root proc 0

	 // Invalidate the input data to save space, since we don't
	 // need it anymore.
	 numEntriesPerRow = null;
	 rowPtr = null;
	 colInd = null;
	 values = null;

	 if (debug && myRank == 0)
	   cerr << "-- Proc 0: About to fill in myRowPtr" << endl;

	 // Allocate and fill in myRowPtr (the row pointer array for
	 // my rank's rows).  We delay this until the end because we
	 // don't need it to compute anything else in distribute().
	 // Each proc can do this work for itself, since it only needs
	 // myNumEntriesPerRow to do so.
	 myRowPtr = arcp<size_type> (myNumRows+1);
	 myRowPtr[0] = 0;
	 for (size_type k = 1; k < myNumRows+1; ++k)
	   myRowPtr[k] = myRowPtr[k-1] + myNumEntriesPerRow[k-1];
	 if (extraDebug && debug)
	   {
	     cerr << "Proc " << Teuchos::rank (*(pRowMap->getComm())) 
		  << ": myRowPtr[0.." << myNumRows << "] = [";
	     for (size_type k = 0; k < myNumRows+1; ++k)
	       {
		 cerr << myRowPtr[k];
		 if (k < myNumRows)
		   cerr << " ";
	       }
	     cerr << "]" << endl << endl;
	   }

	 if (debug && myRank == 0)
	   cerr << "-- Proc 0: Done with distribute" << endl;
      }

      /// \brief Given my proc's data, return the completed sparse matrix.
      /// 
      /// Each proc inserts its data into the sparse matrix, and then,
      /// if callFillComplete is true, all procs call fillComplete().
      /// (For whatever reason, you might not be done with the matrix
      /// yet, so you might want to call fillComplete() yourself.
      /// CrsMatrix::fillResume() doesn't currently work as you might
      /// expect when storage optimization is enabled; it fixes the
      /// graph of the matrix, so that you can't add new entries.)
      /// 
      static sparse_matrix_ptr
      makeMatrix (Teuchos::ArrayRCP<size_t>& myNumEntriesPerRow,
		  Teuchos::ArrayRCP<size_type>& myRowPtr,
		  Teuchos::ArrayRCP<global_ordinal_type>& myColInd,
		  Teuchos::ArrayRCP<scalar_type>& myValues,
		  const map_ptr& pRowMap,
		  const map_ptr& pRangeMap,
		  const map_ptr& pDomainMap,
		  const bool callFillComplete = true)
      {
	using Teuchos::ArrayRCP;
	using Teuchos::ArrayView;
	using Teuchos::rcp;
	using std::cerr;
	using std::endl;

	const bool debug = false;

	// The row pointer array always has at least one entry, even
	// if the matrix has zero rows.  myNumEntriesPerRow, myColInd,
	// and myValues would all be empty arrays in that degenerate
	// case, but the row and domain maps would still be nonnull
	// (though they would be trivial maps).
	TEST_FOR_EXCEPTION(myRowPtr.is_null(), std::logic_error,
			   "makeMatrix: myRowPtr array is null.  "
			   "Please report this bug to the Tpetra developers.");
	TEST_FOR_EXCEPTION(pDomainMap.is_null(), std::logic_error,
			   "makeMatrix: domain map is null.  "
			   "Please report this bug to the Tpetra developers.");
	TEST_FOR_EXCEPTION(pRangeMap.is_null(), std::logic_error,
			   "makeMatrix: range map is null.  "
			   "Please report this bug to the Tpetra developers.");
	TEST_FOR_EXCEPTION(pRowMap.is_null(), std::logic_error,
			   "makeMatrix: row map is null.  "
			   "Please report this bug to the Tpetra developers.");

	// Handy for debugging output; not needed otherwise.
	const int myRank = Teuchos::rank (*(pRangeMap->getComm()));

	if (false && debug)
	  {
	    cerr << "Proc " << myRank << ":" << endl
		 << "-- myRowPtr = [ ";
	    std::copy (myRowPtr.begin(), myRowPtr.end(), 
		       std::ostream_iterator<size_type>(cerr, " "));
	    cerr << "]" << endl << "-- myColInd = [ ";
	    std::copy (myColInd.begin(), myColInd.end(), 
		       std::ostream_iterator<size_type>(cerr, " "));
	    cerr << "]" << endl << endl;
	  }

	// Go through all of my columns, and see if any are not in the
	// domain map.  This is possible if numProcs > 1, otherwise
	// not.
	if (false && debug)
	  {
	    size_type numRemote = 0;
	    std::vector<global_ordinal_type> remoteGIDs;

	    typedef typename ArrayRCP<global_ordinal_type>::const_iterator iter_type;
	    for (iter_type it = myColInd.begin(); it != myColInd.end(); ++it)
	      {
		if (! pDomainMap->isNodeGlobalElement (*it))
		  {
		    numRemote++;
		    remoteGIDs.push_back (*it);
		  }
	      }
	    if (numRemote > 0)
	      {
		cerr << "Proc " << myRank << ": " << numRemote 
		     << " remote GIDs = [ " << endl;
		std::copy (remoteGIDs.begin(), remoteGIDs.end(),
			   std::ostream_iterator<global_ordinal_type>(cerr, " "));
		cerr << "]" << endl;
	      }
	  }

	// Construct the CrsMatrix, using the row map, with the
	// constructor specifying the number of nonzeros for each row.
	// Create with DynamicProfile, so that the fillComplete() can
	// do first-touch reallocation (a NUMA (Non-Uniform Memory
	// Access) optimization on multicore CPUs).
	sparse_matrix_ptr A = 
	  rcp (new sparse_matrix_type (pRowMap, myNumEntriesPerRow, 
				       DynamicProfile));
	TEST_FOR_EXCEPTION(A.is_null(), std::logic_error,
			   "makeMatrix: Initial allocation of CrsMatrix failed"
			   ".  Please report this bug to the Tpetra developers"
			   ".");
	// List of the global indices of my rows.
	// They may or may not be contiguous.
	ArrayView<const global_ordinal_type> myRows = 
	  pRowMap->getNodeElementList();
	const size_type myNumRows = myRows.size();

	// Add this processor's matrix entries to the CrsMatrix.
	for (size_type k = 0; k < myNumRows; ++k)
	  {
	    const size_type myCurPos = myRowPtr[k];
	    const size_t curNumEntries = myNumEntriesPerRow[k];

	    if (false && debug)
	      {
		cerr << "Proc " << myRank << ": k = " << k 
		     << ", myCurPos = " << myCurPos
		     << ", curNumEntries = " << curNumEntries
		     << endl;
	      }
	    // Avoid constructing empty views of ArrayRCP objects.
	    if (curNumEntries > 0)
	      A->insertGlobalValues (myRows[k], 
				     myColInd(myCurPos, curNumEntries),
				     myValues(myCurPos, curNumEntries));
	  }
	// We've entered in all our matrix entries, so we can delete
	// the original data.  This will save memory when we call
	// fillComplete(), so that we never keep more than two copies
	// of the matrix's data in memory at once.
	myNumEntriesPerRow = null;
	myRowPtr = null;
	myColInd = null;
	myValues = null;

	if (callFillComplete)
	  A->fillComplete (pDomainMap, pRangeMap, DoOptimizeStorage);
	return A;
      }

    private:

      /// \brief Read in the Banner line from the given input stream.
      ///
      /// \param in [in/out, valid only on Rank 0] Input stream from 
      ///   which to read the Banner line.
      /// \param lineNumber [in/out, valid only on Rank 0] On input:
      ///   Current line number of the input stream.  On output: if 
      ///   any line(s) were successfully read from the input stream,
      ///   this is incremented by the number of line(s) read.  (This
      ///   includes comment lines.)
      /// \param pComm [in, global] Communicator.
      ///
      /// \return Banner [non-null and valid only on Rank 0]
      static Teuchos::RCP<const Banner>
      readBanner (std::istream& in,
		  size_t& lineNumber,
		  const comm_ptr& pComm,
		  const bool tolerant=false,
		  const bool debug=false)
      {
	using Teuchos::RCP;
	using Teuchos::rcp;
	using std::cerr;
	using std::endl;

	const int myRank = Teuchos::rank (*pComm);

	// The pointer will be non-null on return only on MPI Rank 0.
	// Using a pointer lets the data persist outside the
	// "myRank==0" scopes.
	RCP<Banner> pBanner;
	if (myRank == 0)
	  {
	    typedef Teuchos::ScalarTraits<scalar_type> STS;

	    std::string line;
	    // Try to read a line from the input stream.
	    const bool readFailed = ! getline(in, line);
	    TEST_FOR_EXCEPTION(readFailed, std::invalid_argument,
			       "Failed to get Matrix Market banner line "
			       "from input.");
	    // We read a line from the input stream.
	    lineNumber++; 

	    // Assume that the line we found is the banner line.
	    try {
	      pBanner = rcp (new Banner (line, tolerant));
	    } catch (std::exception& e) {
	      TEST_FOR_EXCEPTION(true, std::invalid_argument, 
				 "Matrix Market banner line contains syntax "
				 "error(s): " << e.what());
	    }
	    // Perform further validation for the special case of
	    // reading in a sparse matrix with entries of type
	    // scalar_type.
	    if (pBanner->matrixType() != "coordinate")
	      throw std::invalid_argument ("Matrix Market input file must contain a "
					   "\"coordinate\"-format sparse matrix in "
					   "order to create a Tpetra::CrsMatrix "
					   "object from it.");
	    else if (! STS::isComplex && pBanner->dataType() == "complex")
	      throw std::invalid_argument ("Matrix Market file contains complex-"
					   "valued data, but your chosen scalar "
					   "type is real.");
	    else if (pBanner->dataType() != "real" && pBanner->dataType() != "complex")
	      throw std::invalid_argument ("Only real or complex data types (no "
					   "pattern or integer matrices) are "
					   "currently supported");
	  }
	return pBanner;
      }

      /// \brief Read sparse matrix dimensions on Rank 0, and broadcast.
      ///
      /// Call on all MPI ranks.  MPI Rank 0 attempts to read in the
      /// coordinate dimensions from the input stream.  If it
      /// succeeds, it broadcasts them to all the other MPI ranks.
      /// (All ranks need to know the matrix dimensions in order to
      /// create domain, range, and column Maps.)
      ///
      /// \return (numRows, numCols, numNonzeros)
      static Teuchos::Tuple<global_ordinal_type, 3>
      readCoordDims (std::istream& in,
		     size_t& lineNumber,
		     const Teuchos::RCP<const Banner>& pBanner,
		     const comm_ptr& pComm,
		     const bool tolerant = false,
		     const bool debug = false)
      {
	// Packed coordinate matrix dimensions (numRows, numCols,
	// numNonzeros); computed on Rank 0 and broadcasted to all MPI
	// ranks.
	Teuchos::Tuple<global_ordinal_type, 3> dims;

	// Read in the coordinate matrix dimensions from the input
	// stream.  "success" tells us whether reading in the
	// coordinate matrix dimensions succeeded ("Guilty unless
	// proven innocent").
	bool success = false;
	if (Teuchos::rank(*pComm) == 0)
	  { 
	    TEST_FOR_EXCEPTION(pBanner->matrixType() != "coordinate", 
			       std::invalid_argument,
			       "The Tpetra::CrsMatrix Matrix Market reader "
			       "only accepts \"coordinate\" (sparse) matrix "
			       "data.");
	    // Unpacked coordinate matrix dimensions
	    global_ordinal_type numRows, numCols, numNonzeros;
	    // Only MPI Rank 0 reads from the input stream
	    success = readCoordinateDimensions (in, numRows, numCols, 
						numNonzeros, lineNumber, 
						tolerant);
	    // Pack up the data into a Tuple so we can send them with
	    // one broadcast instead of three.
	    dims[0] = numRows;
	    dims[1] = numCols;
	    dims[2] = numNonzeros;
	  }
	// Only Rank 0 did the reading, so it decides success.
	//
	// FIXME (mfh 02 Feb 2011) Teuchos::broadcast doesn't know how
	// to send bools.  For now, we convert to/from int instead,
	// using the usual "true is 1, false is 0" encoding.
	{
	  int the_success = success ? 1 : 0; // only matters on MPI Rank 0
	  Teuchos::broadcast (*pComm, 0, &the_success);
	  success = (the_success == 1);
	}
	if (success)
	  // Broadcast (numRows, numCols, numNonzeros) from Rank 0
	  // to all the other MPI ranks.  
	  Teuchos::broadcast (*pComm, 0, dims);
	else
	  // Perhaps in tolerant mode, we could set all the
	  // dimensions to zero for now, and deduce correct
	  // dimensions by reading all of the file's entries and
	  // computing the max(row index) and max(column index).
	  // However, for now we just error out in that case.
	  throw std::invalid_argument ("Error reading Matrix Market sparse "
				       "matrix: failed to read coordinate "
				       "matrix dimensions.");
	return dims;
      }
      
      /// \typedef adder_type
      /// \brief Type of object that "adds" entries to the sparse matrix.
      ///
      /// "Adds" here means that it collects and makes note of matrix
      /// entries read in from the input stream.  This object doesn't
      /// call insertGlobalEntries() or fillComplete() on the
      /// CrsMatrix.  Depending on the Matrix Market banner line
      /// information, it may "symmetrize" the matrix by adding entry
      /// A(j,i) (with the appropriate value depending on the symmetry
      /// type) if entry A(i,j) is seen.
      typedef SymmetrizingAdder<Raw::Adder<scalar_type, global_ordinal_type> > adder_type;

      /// \brief Make an "adder" object for processing matrix data.
      ///
      /// \param pComm [in] Communicator (across whose MPI ranks 
      ///   the sparse matrix will be distributed)
      ///
      /// \param banner [in, nonnull and valid on Rank 0 only]
      ///   Object describing the type and symmetry of matrix data.
      ///
      /// \param dims [in] (numRows, numCols, numEntries).  These are
      ///   the "expected" values as read from the top of the Matrix
      ///   Market input stream.  Whether they are the final values
      ///   depends on the "tolerant" parameter and the actual matrix
      ///   data read from the input stream.
      ///
      /// \param tolerant [in] Whether the adder should be "tolerant"
      ///   of syntax errors and missing/incorrect metadata.  (In
      ///   particular, this refers to the number of rows, columns,
      ///   and entries in the matrix.)
      /// \param debug [in] Whether to print verbose debug output
      ///   to stderr.
      ///
      /// \return An adder_type object [nonnull and valid on Rank 0
      ///   only] that optionally symmetrizes the entries of the
      ///   sparse matrix.
      ///
      static Teuchos::RCP<adder_type>
      makeAdder (const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
		 Teuchos::RCP<const Banner>& pBanner, 
		 const Teuchos::Tuple<global_ordinal_type, 3>& dims,
		 const bool tolerant=false,
		 const bool debug=false)
      {
	if (Teuchos::rank (*pComm) == 0)
	  {
	    using Teuchos::RCP;
	    using Teuchos::rcp;
	    typedef Raw::Adder<scalar_type, global_ordinal_type> raw_adder_type;

	    RCP<raw_adder_type> pRaw (new raw_adder_type (dims[0], dims[1], 
							  dims[2], tolerant, 
							  debug));
	    return rcp (new adder_type (pRaw, pBanner->symmType()));
	  }
	else
	  return Teuchos::null;
      }

    public:

      /// \brief Read sparse matrix from the given Harwell-Boeing file.
      ///
      /// This is a collective operation.  Only Rank 0 opens the file
      /// and reads data from it, but all ranks participate and wait
      /// for the final result.
      /// 
      /// \param filename [in] Name of the Matrix Market file.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param pNode [in] Kokkos Node object.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries 
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static sparse_matrix_ptr
      readSparseFile (const std::string& filename,
		      const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
		      const Teuchos::RCP<node_type>& pNode,
		      const bool callFillComplete=true,
		      const bool tolerant=false,
		      const bool debug=false)
      {
	std::ifstream in (filename.c_str());
	return readSparse (in, pComm, pNode, callFillComplete, tolerant, debug);
      }

      /// \brief Read sparse matrix from the given Harwell-Boeing input stream.
      ///
      /// This is a collective operation.  Only Rank 0 reads data from
      /// the input stream, but all ranks participate and wait for the
      /// final result.
      ///
      /// \param filename [in] The input stream from which to read.
      /// \param pComm [in] Communicator containing all processor(s)
      ///   over which the sparse matrix will be distributed.
      /// \param pNode [in] Kokkos Node object.
      /// \param callFillComplete [in] Whether to call fillComplete()
      ///   on the Tpetra::CrsMatrix, after adding all the entries 
      ///   read in from the input stream.
      /// \param tolerant [in] Whether to read the data tolerantly
      ///   from the file.
      /// \param debug [in] Whether to produce copious status output
      ///   useful for Tpetra developers, but probably not useful for
      ///   anyone else.
      static sparse_matrix_ptr
      readSparse (std::istream& in,	
		  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm, 
		  const Teuchos::RCP<node_type>& pNode,
		  const bool callFillComplete=true,
		  const bool tolerant=false,
		  const bool debug=false)
      {
	using Teuchos::RCP;
	using Teuchos::rcp;
	using Teuchos::Tuple;
	using Teuchos::null;
	using std::cerr;
	using std::endl;
	typedef Teuchos::ScalarTraits<scalar_type> STS;

	const bool extraDebug = false;
	const int myRank = Teuchos::rank (*pComm);

	// Current line number in the input stream.  Various calls
	// will modify this depending on the number of lines that are
	// read from the input stream.
	size_t lineNumber = 1; 	

	if (debug && myRank == 0)
	  cerr << "About to read Matrix Market banner line" << endl;

	// The "Banner" tells you whether the input stream represents
	// a sparse matrix, the symmetry type of the matrix, and the
	// type of the data it contains.
	RCP<const Banner> pBanner = readBanner (in, lineNumber, pComm, 
						tolerant, debug);

	if (debug && myRank == 0)
	  cerr << "About to read Matrix Market dimensions line" << endl;

	// dims = (numRows, numCols, numEntries) (all global),
	// as read from the Matrix Market metadata.
	Tuple<global_ordinal_type, 3> dims = 
	  readCoordDims (in, lineNumber, pBanner, pComm, tolerant, debug);

	if (debug && myRank == 0)
	  cerr << "About to make Adder for collecting matrix data" << endl;

	// "Adder" object for collecting all the sparse matrix entries
	// from the input stream.
	RCP<adder_type> pAdder = 
	  makeAdder (pComm, pBanner, dims, tolerant, debug);

	if (debug && myRank == 0)
	  cerr << "About to read matrix data" << endl;

	// Read the sparse matrix entries from the input stream.
	//
	// "Guilty until proven innocent."
	bool readSuccess = false;
	if (myRank == 0)
	  {
	    // Reader for "coordinate" format sparse matrix data.
	    typedef CoordDataReader<adder_type, global_ordinal_type, 
	      scalar_type, STS::isComplex> reader_type;
	    reader_type reader (pAdder);

	    // Read the sparse matrix entries.
	    //
	    // FIXME (mfh 28 Mar 2011) We should catch exceptions
	    // here, and broadcast any error messages to all
	    // processors so that the exception can be rethrown
	    // everywhere (fail fast and hard, rather than hoping that
	    // MPI propagates the exceptions quickly).
	    std::pair<bool, std::vector<size_t> > results = 
	      reader.read (in, lineNumber, tolerant, debug);

	    readSuccess = results.first;
	  }
	// The broadcast of readSuccess from MPI Rank 0 serves as a
	// barrier.  Note that Tpetra::CrsMatrix::fillComplete() only
	// starts with a barrier in debug mode, so we need some kind
	// of barrier or synchronization beforehand.
	//
	// Currently, Teuchos::broadcast doesn't know how to send
	// bools.  For now, we convert to/from int instead, using the
	// usual "true is 1, false is 0" encoding.
	{
	  int the_readSuccess = readSuccess ? 1 : 0; // only matters on MPI Rank 0
	  Teuchos::broadcast (*pComm, 0, &the_readSuccess);
	  readSuccess = (the_readSuccess == 1);
	}
	// TODO (mfh 01 Feb 2011)
	//
	// In tolerant mode, given a "verbose" flag, report / log any
	// bad line number(s) on MPI Rank 0.
	//
	// Should we run fillComplete() in tolerant mode, if the read
	// did not succeed?
	TEST_FOR_EXCEPTION(! readSuccess, std::runtime_error,
			   "Failed to read in the Matrix Market sparse "
			   "matrix.");
	if (debug && myRank == 0)
	  cerr << "Successfully read the Matrix Market data" << endl;

	// In tolerant mode, we need to rebroadcast the matrix
	// dimensions, since they may be different after reading the
	// actual matrix data.  We only need to broadcast the number
	// of rows and columns.  Only Rank 0 needs to know the actual
	// global number of entries, since (a) we need to merge
	// duplicates on Rank 0 first anyway, and (b) when we
	// distribute the entries, each rank other than Rank 0 will
	// only need to know how many entries it owns, not the total
	// number of entries.
	if (tolerant)
	  {
	    if (debug && myRank == 0)
	      {
		cerr << "Tolerant mode: rebroadcasting matrix dimensions" 
		     << endl
		     << "-- Dimensions before: " << dims[0] << " x " << dims[1]
		     << endl;
	      }
	    // Packed coordinate matrix dimensions (numRows, numCols).
	    Teuchos::Tuple<global_ordinal_type, 2> updatedDims;
	    if (myRank == 0)
	      { // If one or more bottom rows of the matrix contain no
		// entries, then the Adder will report that the number
		// of rows is less than that specified in the
		// metadata.  We allow this case, and favor the
		// metadata so that the zero row(s) will be included.
		updatedDims[0] = std::max (dims[0], pAdder->getAdder()->numRows());
		updatedDims[1] = pAdder->getAdder()->numCols();
	      }
	    Teuchos::broadcast (*pComm, 0, updatedDims);
	    dims[0] = updatedDims[0];
	    dims[1] = updatedDims[1];
	    if (debug && myRank == 0)
	      {
		cerr << "-- Dimensions after: " << dims[0] << " x " << dims[1]
		     << endl;
	      }
	  }
	else 
	  { // In strict mode, we require that the matrix's metadata
	    // and its actual data agree, at least somewhat.  In
	    // particular, the number of rows must agree, since
	    // otherwise we cannot distribute the matrix correctly.

	    // Teuchos::broadcast() doesn't know how to broadcast
	    // bools, so we use an int with the standard 1 == true, 0
	    // == false encoding.
	    int dimsMatch = 1;
	    if (myRank == 0)
	      { // If one or more bottom rows of the matrix contain no
		// entries, then the Adder will report that the number
		// of rows is less than that specified in the
		// metadata.  We allow this case, and favor the
		// metadata, but do not allow the Adder to think there
		// are more rows in the matrix than the metadata says.
		if (dims[0] < pAdder->getAdder()->numRows())
		  dimsMatch = 0;
	      }
	    Teuchos::broadcast (*pComm, 0, &dimsMatch);
	    if (dimsMatch == 0)
	      { // We're in an error state anyway, so we might as well
		// work a little harder to print an informative error
		// message.
		//
		// Broadcast the Adder's idea of the matrix dimensions
		// from Proc 0 to all processes.
		Teuchos::Tuple<global_ordinal_type, 2> addersDims;
		if (myRank == 0)
		  {
		    addersDims[0] = pAdder->getAdder()->numRows();
		    addersDims[1] = pAdder->getAdder()->numCols();
		  }
		Teuchos::broadcast (*pComm, 0, addersDims);
		TEST_FOR_EXCEPTION(dimsMatch == 0, std::runtime_error,
				   "The matrix metadata says that the matrix is "
				   << dims[0] << " x " << dims[1] << ", but the "
				   "actual data says that the matrix is " 
				   << addersDims[0] << " x " << addersDims[1]
				   << ".  That means the data includes more "
				   "rows than reported in the metadata.  This "
				   "is not allowed when parsing in strict mode."
				   "  Parse the matrix in tolerant mode to "
				   "ignore the metadata when it disagrees with "
				   "the data.");
	      }
	  }

	if (debug && myRank == 0)
	  cerr << "Converting matrix data into CSR format on Proc 0" << endl;

	// Now that we've read in all the matrix entries from the
	// input stream into the adder on Rank 0, post-process them
	// into CSR format (still on Rank 0).  This will facilitate
	// distributing them to all the processors.
	//
	// These arrays represent the global matrix data as a CSR
	// matrix (with numEntriesPerRow as redundant but convenient
	// metadata, since it's computable from rowPtr and vice
	// versa).  They are valid only on Rank 0.
	ArrayRCP<size_t> numEntriesPerRow;
	ArrayRCP<size_type> rowPtr;
	ArrayRCP<global_ordinal_type> colInd;
	ArrayRCP<scalar_type> values;

	// Rank 0 first merges duplicate entries, and then converts
	// the coordinate-format matrix data to CSR.
	if (myRank == 0)
	  {
	    typedef Raw::Element<scalar_type, global_ordinal_type> element_type;
	    typedef typename std::vector<element_type>::const_iterator iter_type;

	    if (extraDebug && debug)
	      {
		std::ofstream out ("before.mtx");
		// Print out the read-in matrix entries before the merge.
		const std::vector<element_type>& entries = 
		  pAdder->getAdder()->getEntries();
		std::copy (entries.begin(), entries.end(), 
			   std::ostream_iterator<element_type>(out, "\n"));
	      }
	    // Additively merge duplicate matrix entries.
	    pAdder->getAdder()->merge ();

	    if (extraDebug && debug)
	      {
		std::ofstream out ("after.mtx");
		// Print out the read-in matrix entries after the merge.
		const std::vector<element_type>& entries = 
		  pAdder->getAdder()->getEntries();
		std::copy (entries.begin(), entries.end(), 
			   std::ostream_iterator<element_type>(out, "\n"));
	      }
	    // Get a temporary const view of the merged matrix entries.
	    const std::vector<element_type>& entries = 
	      pAdder->getAdder()->getEntries();

	    // Number of rows in the matrix.  If we are in tolerant
	    // mode, we've already synchronized dims with the actual
	    // matrix data.  If in strict mode, we should use dims (as
	    // read from the file's metadata) rather than the matrix
	    // data to determine the dimensions.  (The matrix data
	    // will claim fewer rows than the metadata, if one or more
	    // rows have no entries stored in the file.)
	    const size_type numRows = dims[0]; 

	    // Number of matrix entries (after merging).
	    const size_type numEntries = entries.size();

	    if (debug)
	      cerr << "Proc 0: Matrix has numRows=" << numRows 
		   << " rows and numEntries=" << numEntries 
		   << " entries." << endl;

	    // Make space for the CSR matrix data.  The conversion to
	    // CSR algorithm is easier if we fill numEntriesPerRow
	    // with zeros at first.
	    numEntriesPerRow = arcp<size_t> (numRows);
	    std::fill (numEntriesPerRow.begin(), numEntriesPerRow.end(), 
		       static_cast<size_t>(0));
	    rowPtr = arcp<size_type> (numRows+1);
	    std::fill (rowPtr.begin(), rowPtr.end(), 
		       static_cast<size_type>(0));
	    colInd = arcp<global_ordinal_type> (numEntries);
	    values = arcp<scalar_type> (numEntries);

	    // Convert from array-of-structs coordinate format to CSR
	    // (compressed sparse row) format.
	    {
	      global_ordinal_type prvRow = 0;
	      size_type curPos = 0;
	      rowPtr[0] = 0;
	      for (curPos = 0; curPos < numEntries; ++curPos)
		{
		  const element_type& curEntry = entries[curPos];
		  const global_ordinal_type curRow = curEntry.rowIndex();
		  TEST_FOR_EXCEPTION(curRow < prvRow, std::logic_error,
				     "Row indices out of order, even though "
				     "they are supposed to be sorted.  curRow"
				     " = " << curRow << ", prvRow = " 
				     << prvRow << ", at curPos = " << curPos 
				     << ".  Please report this bug to the "
				     "Tpetra developers.");
		  if (curRow > prvRow)
		    {
		      for (size_type r = prvRow+1; r <= curRow; ++r)
			rowPtr[r] = curPos;
		      prvRow = curRow;
		    }
		  numEntriesPerRow[curRow]++;
		  colInd[curPos] = curEntry.colIndex();
		  values[curPos] = curEntry.value();
		}
	      // rowPtr has one more entry than numEntriesPerRow.  The
	      // last entry of rowPtr is the number of entries in colInd
	      // and values.
	      rowPtr[numRows] = numEntries;
	    }
	    if (debug)
	      {
		const size_type maxToDisplay = 100;

		cerr << "Proc 0: numEntriesPerRow[0.." << (numEntriesPerRow.size()-1) << "] ";
		if (numRows > maxToDisplay)
		  cerr << "(only showing first and last few entries) ";
		cerr << "= [";
		if (numRows > 0)
		  {
		    if (numRows > maxToDisplay)
		      {
			for (size_type k = 0; k < 2; ++k)
			  cerr << numEntriesPerRow[k] << " ";
			cerr << "... ";
			for (size_type k = numRows-2; k < numRows-1; ++k)
			  cerr << numEntriesPerRow[k] << " ";
		      }
		    else
		      {		   
			for (size_type k = 0; k < numRows-1; ++k)
			  cerr << numEntriesPerRow[k] << " ";
		      }
		    cerr << numEntriesPerRow[numRows-1];
		  }
		cerr << "]" << endl;

		cerr << "Proc 0: rowPtr ";
		if (numRows > maxToDisplay)
		  cerr << "(only showing first and last few entries) ";
		cerr << "= [";
		if (numRows > maxToDisplay)
		  {
		    for (size_type k = 0; k < 2; ++k)
		      cerr << rowPtr[k] << " ";
		    cerr << "... ";
		    for (size_type k = numRows-2; k < numRows; ++k)
		      cerr << rowPtr[k] << " ";
		  }
		else
		  {
		    for (size_type k = 0; k < numRows; ++k)
		      cerr << rowPtr[k] << " ";
		  }
		cerr << rowPtr[numRows] << "]" << endl;
	      }
	  }
	// Now we're done with the Adder, so we can release the
	// reference ("free" it) to save space.
	pAdder = null;

	if (debug && myRank == 0)
	  cerr << "Making range, domain, and row maps" << endl;

	// Make the maps that describe the matrix'x range and domain,
	// and the distribution of its rows.
	map_ptr pRangeMap = makeRangeMap (pComm, pNode, dims[0]);
	map_ptr pDomainMap = makeDomainMap (pRangeMap, dims[0], dims[1]);
	map_ptr pRowMap = makeRowMap (null, pComm, pNode, dims[0]);

	if (debug && myRank == 0)
	  cerr << "Distributing the matrix data" << endl;

	// Distribute the matrix data.  Each processor has to add the
	// rows that it owns.  If you try to make Rank 0 call
	// insertGlobalValues() for _all_ the rows, not just those it
	// owns, then fillComplete() will compute the number of
	// columns incorrectly.  That's why Rank 0 has to distribute
	// the matrix data and why we make all the processors (not
	// just Rank 0) call insertGlobalValues() on their own data.
	//
	// These arrays represent each processor's part of the matrix
	// data, in "CSR" format (sort of, since the row indices might
	// not be contiguous).
	ArrayRCP<size_t> myNumEntriesPerRow;
	ArrayRCP<size_type> myRowPtr;
	ArrayRCP<global_ordinal_type> myColInd;
	ArrayRCP<scalar_type> myValues;
	// Distribute the matrix data.
	distribute (myNumEntriesPerRow, myRowPtr, myColInd, myValues, pRowMap,
		    numEntriesPerRow, rowPtr, colInd, values, debug);

	if (debug && myRank == 0)
	  {
	    cerr << "Inserting matrix entries on each processor";
	    if (callFillComplete)
	      cerr << " and calling fillComplete()";
	    cerr << endl;
	  }
	// Each processor inserts its part of the matrix data, and
	// then they all call fillComplete().  This method invalidates
	// the my* distributed matrix data before calling
	// fillComplete(), in order to save space.  In general, we
	// never store more than two copies of the matrix's entries in
	// memory at once, which is no worse than what Tpetra
	// promises.
	sparse_matrix_ptr pMatrix = 
	  makeMatrix (myNumEntriesPerRow, myRowPtr, myColInd, myValues,
		      pRowMap, pRangeMap, pDomainMap, callFillComplete);
	TEST_FOR_EXCEPTION(pMatrix.is_null(), std::logic_error,
			   "makeMatrix() returned a null pointer.  Please "
			   "report this bug to the Tpetra developers.");

	// We can't get the dimensions of the matix until after
	// fillComplete() is called.  Thus, we can't do the sanity
	// check (dimensions read from the Matrix Market data,
	// vs. dimensions reported by the CrsMatrix) unless the user
	// asked makeMatrix() to call fillComplete().
	//
	// Note that pMatrix->getGlobalNum{Rows,Cols}() does _not_ do
	// what one might think it does, so you have to ask the range
	// resp. domain map for the number of rows resp. columns.
	if (callFillComplete)
	  {
	    if (extraDebug && debug)
	      {
		const global_size_t globalNumRows = 
		  pRangeMap->getGlobalNumElements();
		const global_size_t globalNumCols = 
		  pDomainMap->getGlobalNumElements();
		if (myRank == 0)
		  {
		    cerr << "-- Matrix is " 
			 << globalNumRows << " x " << globalNumCols 
			 << " with " << pMatrix->getGlobalNumEntries()
			 << " entries, and index base " 
			 << pMatrix->getIndexBase() << "." << endl;
		  }
		Teuchos::barrier (*pComm);
		for (int p = 0; p < Teuchos::size (*pComm); ++p)
		  {
		    if (myRank == p)
		      {
			cerr << "-- Proc " << p << " owns " 
			     << pMatrix->getNodeNumCols() 
			     << " columns, and " 
			     << pMatrix->getNodeNumEntries()
			     << " entries." << endl;
		      }
		    Teuchos::barrier (*pComm);
		  }
	      } // if (extraDebug && debug)
	  } // if (callFillComplete)

	if (debug && myRank == 0)
	  cerr << "Done creating the Tpetra::CrsMatrix from the Matrix Market "
	    "data" << endl;

	return pMatrix;
      }
    };



    /// \class Writer
    /// \brief Matrix Market file writer for Tpetra::CrsMatrix.
    ///
    /// The writeSparse() and writeSparseFile() class methods write a
    /// Tpetra::CrsMatrix sparse matrix to a Matrix Market
    /// "coordinate" format sparse matrix output stream resp. file.
    /// Currently, they assume that all ranks can write to the output
    /// stream; this may be fixed in the future, so that they only
    /// write to the output stream on Rank 0.
    template<class SparseMatrixType>
    class Writer {
    public:
      typedef SparseMatrixType sparse_matrix_type;
      typedef Teuchos::RCP<sparse_matrix_type> sparse_matrix_ptr;
      /// \typedef scalar_type
      /// \brief Type of the entries of the sparse matrix.
      typedef typename SparseMatrixType::scalar_type scalar_type;
      /// \typedef local_ordinal_type
      /// \brief Only used to define map_type.
      typedef typename SparseMatrixType::local_ordinal_type local_ordinal_type;
      /// \typedef global_ordinal_type
      /// \brief Type of indices as read from the Matrix Market file.
      ///
      /// Indices of the sparse matrix are stored as global ordinals,
      /// since Matrix Market files represent the whole matrix and
      /// don't have a notion of distribution.
      typedef typename SparseMatrixType::global_ordinal_type global_ordinal_type;
      typedef typename SparseMatrixType::node_type node_type;

      static void
      writeSparseFile (const std::string& filename,
		       Teuchos::RCP<const sparse_matrix_type>& pMatrix)
      {
	TEST_FOR_EXCEPTION(true, std::logic_error, 
			   "Not implemented yet.  First, we need to change "
			   "writeSparse() so that it only writes to the output "
			   "stream argument on Rank 0.");
      }

      /// \brief Print the sparse matrix in Matrix Market format.
      ///
      /// Write the given Tpetra::CrsMatrix sparse matrix to a Matrix
      /// Market "coordinate" format sparse matrix output stream.
      /// Currently, this method assumes that all ranks can write to
      /// the output stream; this will be fixed in the future, so that
      /// it only writes to the output stream on Rank 0.
      static void
      writeSparse (std::ostream& out,
		   const Teuchos::RCP<const sparse_matrix_type>& pMatrix)
      {
	using Teuchos::Comm;
	using Teuchos::RCP;
	using std::endl;
	typedef typename Teuchos::ScalarTraits<scalar_type> STS;
	typedef typename STS::magnitudeType magnitude_type;
	typedef typename Teuchos::ScalarTraits<magnitude_type> STM;
	typedef typename ArrayView<scalar_type>::size_type size_type;

	RCP<const Comm<int> > pComm = pMatrix->getComm();
	const int numProcs = Teuchos::size (*pComm);
	const int myRank = Teuchos::rank (*pComm);

	//
	// Set up the output stream.
	//
	// Print floating-point values in scientific notation.
	out << std::scientific;

	// We're writing decimal digits, so compute the number of
	// digits we need to get reasonable accuracy when reading
	// values back in.
	//
	// There is actually an algorithm, due to Guy Steele (yes,
	// _that_ Guy Steele) et al., for idempotent printing of
	// finite-length floating-point values.  We should actually
	// implement that algorithm, but I don't have time for that
	// now.  Currently, I just print no more than (one decimal
	// digit more than (the number of decimal digits justified by
	// the precision of magnitude_type)).
	{
	  // FIXME (mfh 06 Apr 2011) This will only work if log10 is
	  // defined for magnitude_type inputs.  Teuchos::ScalarTraits
	  // does not currently have an log10() class method.
	  const magnitude_type numDecDigits = 
	    STM::t() * std::log10 (STM::base());

	  // Round and add one. Hopefully this doesn't overflow...
	  const magnitude_type one = STM::one();
	  const magnitude_type two = one + one;
	  const int prec = 1 + static_cast<int> ((two*numDecDigits + one) / two);

	  // Set the number of (decimal) digits after the decimal
	  // point to print.
	  out.precision (prec);
	}	

	// Number of errors encountered while printing out the matrix.
	int numPrintErrors = 0;

	// Rank 0: Print the banner line and the dimensions line to out.
	if (myRank == 0)
	  {
	    out << "%%MatrixMarket matrix coordinate ";
	    if (STS::isComplex)
	      out << "complex ";
	    else
	      out << "real ";
	    out << "general" << endl;

	    // getGlobalNum{Rows,Cols}() does not return what you
	    // think it should return.  Instead, ask the range
	    // resp. domain map for the number of rows
	    // resp. columns in the matrix.
	    out << pMatrix->getRangeMap()->getGlobalNumElements() 
		<< " "
		<< pMatrix->getDomainMap()->getGlobalNumElements() 
		<< " "
		<< pMatrix->getGlobalNumEntries()
		<< endl;
	  }
	Teuchos::barrier (*pComm);

	// Let each processor in turn print its rows to out.  We
	// assume that all processors can print to out.  We do
	// _not_ assume here that the row map is one-to-one;
	// printing should work just fine, as long as nonzeros
	// themselves are not stored redundantly.
	for (int p = 0; p < numProcs; ++p)
	  {
	    if (myRank == p)
	      {
		// Storage for column indices and values in each
		// row.  Will be resized as necessary.  (Proc p
		// may not own any rows, in which case Proc p
		// won't need to allocate these at all.
		Array<global_ordinal_type> indices;
		Array<scalar_type> values;
		
		// List of the rows with storage on Proc p.
		ArrayView<const global_ordinal_type> myRows = 
		  pMatrix->getRowMap()->getNodeElementList();
		// Number of rows with storage on Proc p.
		const size_type myNumRows = myRows.size();
		
		// For each row that Proc p owns, print its entries.
		for (size_type k = 0; k < myNumRows; ++k)
		  {
		    const global_ordinal_type curRow = myRows[k];
		    size_t numEntries = 
		      pMatrix->getNumEntriesInGlobalRow (curRow);
		    
		    // Resize (if necessary) the arrays for holding
		    // column indices and values for the current row.
		    //
		    // Signed to unsigned integer conversion, for
		    // integers of the same size, shouldn't overflow.
		    if (static_cast<size_t> (indices.size()) < numEntries)
		      indices.resize (numEntries);
		    if (static_cast<size_t> (values.size()) < numEntries)
		      values.resize (numEntries);
		    // This views are exactly the right length to hold
		    // the data for the current row.  indices and
		    // values may be longer than necessary; that's an
		    // optimization, to avoid resizing them with every
		    // row.
		    ArrayView<global_ordinal_type> indicesView = 
		      indices.view (0, numEntries);
		    ArrayView<scalar_type> valuesView = 
		      values.view (0, numEntries);

		    // Make sure there were no surprises with the
		    // number of entries.
		    size_t newNumEntries = 0;
		    pMatrix->getGlobalRowCopy (curRow, indicesView,
					       valuesView, newNumEntries);
		    if (newNumEntries != numEntries)
		      numPrintErrors++;
		    else
		      {
			// This will convert the current row index to
			// the one-based indices that Matrix Market
			// files want.  Ditto for the column index
			// code within the for loop.
			const global_ordinal_type indexBase = 
			  pMatrix->getIndexBase();
			const global_ordinal_type curRowIdx = 
			  curRow +
			  static_cast<global_ordinal_type>(1) - 
			  indexBase;
			// Loop over all entries in the current row,
			// printing each entry on a separate line.
			for (size_t j = 0; j < numEntries; ++j)
			  {
			    const global_ordinal_type curColIdx = 
			      indicesView[j] + 
			      static_cast<global_ordinal_type>(1) - 
			      indexBase;
			    out << curRowIdx << " " 
				<< curColIdx << " ";
			    if (STS::isComplex)
			      out << STS::real(valuesView[j]) << " " 
				  << STS::imag(valuesView[j]);
			    else
			      out << valuesView[j];
			    out << endl;
			  }
		      }
		  }
		Teuchos::barrier (*pComm);

		// If there were any errors on any processors, stop
		// right away and throw an exception.
		int totalNumPrintErrors = 0;
		Teuchos::reduceAll (*pComm, Teuchos::REDUCE_SUM, numPrintErrors,
				    Teuchos::Ptr<int> (&totalNumPrintErrors));
		TEST_FOR_EXCEPTION(totalNumPrintErrors > 0, std::runtime_error,
				   "Failed to print Tpetra::CrsMatrix.  Total "
				   "number of print errors thus far: " 
				   << numPrintErrors);
	      } // If myRank == p
	  } // For each proc p in 0, 1, ..., numProcs-1
      }
    }; // class Writer
    
  } // namespace MatrixMarket
} // namespace Tpetra

#endif // __MatrixMarket_Tpetra_hpp
