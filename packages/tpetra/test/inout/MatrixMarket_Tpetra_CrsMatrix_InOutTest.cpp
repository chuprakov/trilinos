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

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Util.hpp> // sort2
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_GlobalMPISession.hpp>

using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::as;
using Teuchos::Comm;
using Teuchos::OSTab;
using Teuchos::ParameterList;
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MAX;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;

const bool callFillComplete = true;
const bool tolerant = false;
// Whether to print copious debugging output to stderr when doing
// Matrix Market input and output. 
const bool debug = false;

namespace {

const char matrix_symRealSmall[] = 
"%%MatrixMarket matrix coordinate real general\n"
"5 5 13\n"
"1 1  2.0000000000000e+00\n"
"1 2  -1.0000000000000e+00\n"
"2 1  -1.0000000000000e+00\n"
"2 2  2.0000000000000e+00\n"
"2 3  -1.0000000000000e+00\n"
"3 2  -1.0000000000000e+00\n"
"3 3  2.0000000000000e+00\n"
"3 4  -1.0000000000000e+00\n"
"4 3  -1.0000000000000e+00\n"
"4 4  2.0000000000000e+00\n"
"4 5  -1.0000000000000e+00\n"
"5 4  -1.0000000000000e+00\n"
"5 5  2.0000000000000e+00\n";

// Input matrices must be fill complete.
template<class CrsMatrixType>
bool
compareCrsMatrix (const CrsMatrixType& A_orig, const CrsMatrixType& A)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  typedef typename CrsMatrixType::scalar_type ST;
  typedef typename CrsMatrixType::global_ordinal_type GO;

  if (! A_orig.getRowMap ()->isSameAs (* (A.getRowMap ()))) {
    return false;
  }
  else if (! A_orig.getColMap ()->isSameAs (* (A.getColMap ()))) {
    return false;
  }
  else if (! A_orig.getDomainMap ()->isSameAs (* (A.getDomainMap ()))) {
    return false;
  }
  else if (! A_orig.getRangeMap ()->isSameAs (* (A.getRangeMap ()))) {
    return false;
  }
  else {
    //
    // Are my local matrices equal?
    //
    RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
    int localEqual = 1;

    Array<GO> indOrig, ind;
    Array<ST> valOrig, val;
    size_t numEntriesOrig = 0;
    size_t numEntries = 0;

    ArrayView<const GO> localElts = A.getRowMap ()->getNodeElementList ();
    typedef typename ArrayView<const GO>::size_type size_type;
    const size_type numLocalElts = localElts.size ();
    for (size_type k = 0; k < numLocalElts; ++k) {
      const GO globalRow = localElts[k];
      numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
      numEntries = A.getNumEntriesInGlobalRow (globalRow);

      if (numEntriesOrig != numEntries) {
	localEqual = 0;
	break;
      }
      indOrig.resize (numEntriesOrig);
      valOrig.resize (numEntriesOrig);      
      A_orig.getGlobalRowCopy (globalRow, indOrig (), valOrig (), numEntriesOrig);
      ind.resize (numEntries);
      val.resize (numEntries);      
      A.getGlobalRowCopy (globalRow, ind (), val (), numEntries);

      // Global row entries are not necessarily sorted.  Sort them so
      // we can compare them.
      Tpetra::sort2 (indOrig.begin (), indOrig.end (), valOrig.begin ());
      Tpetra::sort2 (ind.begin (), ind.end (), val.begin ());

      for (size_t entryIndex = 0; entryIndex < numEntries; ++entryIndex) {
	// Values should be _exactly_ equal.
	if (indOrig[k] != ind[k] || valOrig[k] != val[k]) {
	  localEqual = 0;
	  break;
	}
      }
    }

    int globalEqual = 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
    return globalEqual == 1;
  }
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> >
createSymRealSmall (const Teuchos::RCP<const Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, NodeType> >& rowMap)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm ();
  const int myRank = comm->getRank ();

  const GST globalNumElts = rowMap->getGlobalNumElements ();
  const size_t myNumElts = (myRank == 0) ? as<size_t> (globalNumElts) : 0;
  const GO indexBase = rowMap->getIndexBase ();
  RCP<const map_type> gatherRowMap = 
    rcp (new map_type (globalNumElts, myNumElts, indexBase, 
		       comm, rowMap->getNode ()));
  matrix_type A_gather (gatherRowMap, as<size_t> (0));

  if (myRank == 0) {
    Array<GO> ind (3);
    Array<ST> val (3);
    for (size_t myRow = 0; myRow < myNumElts; ++myRow) {
      const GO globalRow = gatherRowMap->getGlobalElement (myRow);
      if (globalRow == gatherRowMap->getMinAllGlobalIndex ()) {
	val[0] = as<ST> (2);
	val[1] = as<ST> (-1);
	ind[0] = globalRow;
	ind[1] = globalRow + 1;
	A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else if (globalRow == gatherRowMap->getMaxAllGlobalIndex ()) {
	val[0] = as<ST> (-1);
	val[1] = as<ST> (2);
	ind[0] = globalRow - 1;
	ind[1] = globalRow;
	A_gather.insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
      }
      else {
	val[0] = as<ST> (-1);
	val[1] = as<ST> (2);
	val[2] = as<ST> (-1);
	ind[0] = globalRow - 1;
	ind[1] = globalRow;
	ind[2] = globalRow + 1;
	A_gather.insertGlobalValues (globalRow, ind.view (0, 3), val.view (0, 3));
      }
    }
  }
  A_gather.fillComplete (rowMap, rowMap);
  RCP<matrix_type> A = rcp (new matrix_type (rowMap, as<size_t> (0)));
  export_type exp (gatherRowMap, rowMap);
  A->doExport (A_gather, exp, Tpetra::INSERT);
  A->fillComplete (rowMap, rowMap);
  return A;
}

template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
void
testCrsMatrix (Teuchos::FancyOStream& out, const GlobalOrdinalType indexBase)
{
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> crs_matrix_type;
  bool result = true; // current Boolean result; reused below
  bool success = true; // used by TEST_EQUALITY

  out << "Test: CrsMatrix Matrix Market I/O, w/ Map with index base " 
      << indexBase << endl;
  OSTab tab1 (out);

  RCP<const Comm<int> > comm = 
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();

  out << "Original sparse matrix:" << endl;
  out << matrix_symRealSmall << endl;

  out << "Creating the row Map" << endl;
  const global_size_t globalNumElts = 5;
  RCP<const map_type> rowMap = 
    rcp (new map_type (globalNumElts, indexBase, comm, 
		       Tpetra::GloballyDistributed, node));

  out << "Reading in the matrix" << endl;
  std::istringstream inStr (matrix_symRealSmall);
  RCP<const map_type> colMap;
  RCP<const map_type> domainMap = rowMap;
  RCP<const map_type> rangeMap = rowMap;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  RCP<crs_matrix_type> A = 
    reader_type::readSparse (inStr, rowMap, colMap, domainMap, rangeMap, 
			     callFillComplete, tolerant, debug);

  out << "Creating original matrix" << endl;
  RCP<crs_matrix_type> A_orig = createSymRealSmall<ST, LO, GO, NT> (rowMap);

  out << "Comparing read-in matrix to original matrix" << endl;
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A);
  TEST_EQUALITY( result, true );

  out << "Writing out the original matrix" << endl;
  std::ostringstream outStr;
  typedef Tpetra::MatrixMarket::Writer<crs_matrix_type> writer_type;
  writer_type::writeSparse (outStr, A_orig, debug);

  out << "Reading it in again and comparing with original matrix" << endl;
  std::istringstream inStr2 (outStr.str ());
  RCP<crs_matrix_type> A_orig2 = 
    reader_type::readSparse (inStr2, rowMap, colMap, domainMap, rangeMap, 
			     callFillComplete, tolerant, debug);
  result = compareCrsMatrix<crs_matrix_type> (*A_orig, *A_orig2);
  TEST_EQUALITY( result, true );
}

} // namespace (anonymous)

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase0, ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  const GlobalOrdinalType indexBase = 0;
  testCrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> (out, indexBase);
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrixOutputInput, IndexBase1, ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType )
{
  const GlobalOrdinalType indexBase = 1;
  testCrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> (out, indexBase);
}


// Unit test macro isn't smart enough to deal with namespace qualifications.
typedef Kokkos::DefaultNode::DefaultNodeType the_node_type;

// We instantiate tests for all combinations of the following parameters:
// - indexBase = {0, 1}
// - ST = {double, float}
// - GO = {int, long}
//
// We should really use the Tpetra ETI system to control which GO we
// test here, but int and long are the two most important cases.

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, double, int, int, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, double, int, int, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase0, double, int, long, the_node_type )

TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrixOutputInput, IndexBase1, double, int, long, the_node_type )

