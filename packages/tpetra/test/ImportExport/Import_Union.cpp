/*
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
*/

#include "Teuchos_UnitTestHarness.hpp"

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Distributor.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Map.hpp>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>

#include <algorithm>


namespace {

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::global_size_t GST;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions (true);
  }

  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( ImportUnion, ContigPlusContig, LocalOrdinalType, GlobalOrdinalType )
  {
    //typedef double ST;
    typedef LocalOrdinalType LO;
    typedef GlobalOrdinalType GO;
    typedef Kokkos::SerialNode NT;
    typedef Tpetra::Map<LO, GO, NT> map_type;
    typedef Tpetra::Import<LO, GO, NT> import_type;
    typedef typename Array<GO>::size_type size_type;

    out << "Tpetra::Import::setUnion test" << endl;
    OSTab tab1 (out);
    out << "Both target Maps contiguous" << endl;

    RCP<const Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
    RCP<NT> node = Tpetra::DefaultPlatform::getDefaultPlatform ().getNode ();
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const GO indexBase = 0;
    const int myRank = comm->getRank ();
    const int numProcs = comm->getSize ();

    const GO n1 = 5;
    const GO n2 = 10;

    // Let r = comm->getRank() and P = comm->getSize().  Then:
    //
    // Target Map 1: indexBase + {min(n1*r - 1, 0), ..., max(n1*r + n1, n1*P - 1)}.
    // Target Map 2: indexBase + {min(n2*r - 2, 0), ..., max(n2*r + n2, n2*P - 2)}.
    // Source Map (for both): indexBase + {0, ..., max(n1*r + n1 - 1, n1*P - 1)}.

    Array<GO> srcMapGids;
    for (GO k = 0; k < std::max (n1*myRank + n1, n1*numProcs); ++k) {
      srcMapGids.push_back (indexBase + k);
    }
    Array<GO> tgtMap1Gids;
    // WARNING (mfh 24 Apr 2013) This ONLY works if GO is signed.
    for (GO k = std::max (as<GO> (n1*myRank - 1), as<GO> (0));
         k < std::min (as<GO> (n1*myRank + n1 + 1), as<GO> (n1*numProcs));
         ++k) {
      tgtMap1Gids.push_back (indexBase + k);
    }
    Array<GO> tgtMap2Gids;
    // WARNING (mfh 24 Apr 2013) This ONLY works if GO is signed.
    for (size_type k = std::max (as<GO> (n2*myRank - 1), as<GO> (0));
         k < std::min (as<GO> (n2*myRank + n2 + 1), as<GO> (n2*numProcs));
         ++k) {
      tgtMap2Gids.push_back (indexBase + k);
    }
    Array<GO> unionMapGids;
    Array<GO> tgtMap1GidsCopy (tgtMap1Gids.begin (), tgtMap1Gids.end ());
    Array<GO> tgtMap2GidsCopy (tgtMap2Gids.begin (), tgtMap2Gids.end ());

    std::sort (tgtMap1GidsCopy.begin (), tgtMap1GidsCopy.end ());
    std::sort (tgtMap2GidsCopy.begin (), tgtMap2GidsCopy.end ());
    std::set_union (tgtMap1GidsCopy.begin (), tgtMap1GidsCopy.end (),
                    tgtMap2GidsCopy.begin (), tgtMap2GidsCopy.end (),
                    std::back_inserter (unionMapGids));

    out << "Making the Maps" << endl;

    RCP<const map_type> srcMap (new map_type (INVALID, srcMapGids (), indexBase, comm, node));
    RCP<const map_type> tgtMap1 (new map_type (INVALID, tgtMap1GidsCopy (), indexBase, comm, node));
    RCP<const map_type> tgtMap2 (new map_type (INVALID, tgtMap2GidsCopy (), indexBase, comm, node));
    RCP<const map_type> expectedUnionMap (new map_type (INVALID, unionMapGids (), indexBase, comm, node));

    out << "Making the Import objects" << endl;

    RCP<const import_type> imp1 (new import_type (srcMap, tgtMap1));
    RCP<const import_type> imp2 (new import_type (srcMap, tgtMap2));
    RCP<const import_type> expectedUnionImp (new import_type (srcMap, expectedUnionMap));

    out << "Computing setUnion using first, second" << endl;
    RCP<const import_type> unionImp1 = imp1->setUnion (*imp2);

    out << "Computing setUnion using second, first" << endl;
    RCP<const import_type> unionImp2 = imp2->setUnion (*imp1);

    out << "Running tests" << endl;

    const bool targetMapSame1 = expectedUnionMap->isSameAs (* (unionImp1->getTargetMap ()));
    TEST_EQUALITY( targetMapSame1, true );

    const bool targetMapSame2 = expectedUnionMap->isSameAs (* (unionImp2->getTargetMap ()));
    TEST_EQUALITY( targetMapSame2, true );

    const bool numSameIDsSame1 = unionImp1->getNumSameIDs () == expectedUnionImp->getNumSameIDs ();
    TEST_EQUALITY( numSameIDsSame1, true );

    const bool numSameIDsSame2 = unionImp2->getNumSameIDs () == expectedUnionImp->getNumSameIDs ();
    TEST_EQUALITY( numSameIDsSame2, true );

    const bool numPermuteIDsSame1 = unionImp1->getNumPermuteIDs () == expectedUnionImp->getNumPermuteIDs ();
    TEST_EQUALITY( numPermuteIDsSame1, true );

    const bool numPermuteIDsSame2 = unionImp2->getNumPermuteIDs () == expectedUnionImp->getNumPermuteIDs ();
    TEST_EQUALITY( numPermuteIDsSame2, true );

    const bool numRemoteIDsSame1 = unionImp1->getNumRemoteIDs () == expectedUnionImp->getNumRemoteIDs ();
    TEST_EQUALITY( numRemoteIDsSame1, true );

    const bool numRemoteIDsSame2 = unionImp2->getNumRemoteIDs () == expectedUnionImp->getNumRemoteIDs ();
    TEST_EQUALITY( numRemoteIDsSame2, true );

    const bool numExportIDsSame1 = unionImp1->getNumExportIDs () == expectedUnionImp->getNumExportIDs ();
    TEST_EQUALITY( numExportIDsSame1, true );

    const bool numExportIDsSame2 = unionImp2->getNumExportIDs () == expectedUnionImp->getNumExportIDs ();
    TEST_EQUALITY( numExportIDsSame2, true );

    ArrayView<const LO> permuteFromLIDs_actual1 = unionImp1->getPermuteFromLIDs ();
    ArrayView<const LO> permuteFromLIDs_actual2 = unionImp2->getPermuteFromLIDs ();
    ArrayView<const LO> permuteFromLIDs_expected = expectedUnionImp->getPermuteFromLIDs ();
    const bool permuteFromLIDsSame1 =
      numPermuteIDsSame1 &&
      permuteFromLIDs_actual1.size () == permuteFromLIDs_expected.size () &&
      std::equal (permuteFromLIDs_actual1.begin (),
                  permuteFromLIDs_actual1.end (),
                  permuteFromLIDs_expected.begin ());
    TEST_EQUALITY( permuteFromLIDsSame1, true );
    const bool permuteFromLIDsSame2 =
      numPermuteIDsSame2 &&
      permuteFromLIDs_actual2.size () == permuteFromLIDs_expected.size () &&
      std::equal (permuteFromLIDs_actual2.begin (),
                  permuteFromLIDs_actual2.end (),
                  permuteFromLIDs_expected.begin ());
    TEST_EQUALITY( permuteFromLIDsSame2, true );

    ArrayView<const LO> permuteToLIDs_actual1 = unionImp1->getPermuteToLIDs ();
    ArrayView<const LO> permuteToLIDs_actual2 = unionImp2->getPermuteToLIDs ();
    ArrayView<const LO> permuteToLIDs_expected = expectedUnionImp->getPermuteToLIDs ();
    const bool permuteToLIDsSame1 =
      numPermuteIDsSame1 &&
      permuteToLIDs_actual1.size () == permuteToLIDs_expected.size () &&
      std::equal (permuteToLIDs_actual1.begin (),
                  permuteToLIDs_actual1.end (),
                  permuteToLIDs_expected.begin ());
    TEST_EQUALITY( permuteToLIDsSame1, true );
    const bool permuteToLIDsSame2 =
      numPermuteIDsSame2 &&
      permuteToLIDs_actual2.size () == permuteToLIDs_expected.size () &&
      std::equal (permuteToLIDs_actual2.begin (),
                  permuteToLIDs_actual2.end (),
                  permuteToLIDs_expected.begin ());
    TEST_EQUALITY( permuteToLIDsSame2, true );

    Tpetra::Distributor& dist_expected = expectedUnionImp->getDistributor ();
    Tpetra::Distributor& dist_actual1 = unionImp1->getDistributor ();
    Tpetra::Distributor& dist_actual2 = unionImp2->getDistributor ();

    const bool sameNumReceives1 =
      dist_expected.getNumReceives () == dist_actual1.getNumReceives ();
    TEST_EQUALITY( sameNumReceives1, true );
    const bool sameNumSends1 =
      dist_expected.getNumSends () == dist_actual1.getNumSends ();
    TEST_EQUALITY( sameNumSends1, true );
    const bool sameHasSelfMessage1 =
      dist_expected.hasSelfMessage () == dist_actual1.hasSelfMessage ();
    TEST_EQUALITY( sameHasSelfMessage1, true );
    const bool sameMaxSendLength1 =
      dist_expected.getMaxSendLength () == dist_actual1.getMaxSendLength ();
    TEST_EQUALITY( sameMaxSendLength1, true );
    const bool sameTotalReceiveLength1 =
      dist_expected.getTotalReceiveLength () == dist_actual1.getTotalReceiveLength ();
    TEST_EQUALITY( sameTotalReceiveLength1, true );

    ArrayView<const int> imagesFrom_expected = dist_expected.getImagesFrom ();
    ArrayView<const int> imagesFrom_actual1 = dist_actual1.getImagesFrom ();
    const bool sameImagesFrom1 =
      imagesFrom_expected.size () == imagesFrom_actual1.size () &&
      std::equal (imagesFrom_expected.begin (),
                  imagesFrom_expected.end (),
                  imagesFrom_actual1.begin ());
    TEST_EQUALITY( sameImagesFrom1, true );

    ArrayView<const int> imagesTo_expected = dist_expected.getImagesTo ();
    ArrayView<const int> imagesTo_actual1 = dist_actual1.getImagesTo ();
    const bool sameImagesTo1 =
      imagesTo_expected.size () == imagesTo_actual1.size () &&
      std::equal (imagesTo_expected.begin (),
                  imagesTo_expected.end (),
                  imagesTo_actual1.begin ());
    TEST_EQUALITY( sameImagesTo1, true );

    ArrayView<const size_t> lengthsFrom_expected = dist_expected.getLengthsFrom ();
    ArrayView<const size_t> lengthsFrom_actual1 = dist_actual1.getLengthsFrom ();
    const bool sameLengthsFrom1 =
      lengthsFrom_expected.size () == lengthsFrom_actual1.size () &&
      std::equal (lengthsFrom_expected.begin (),
                  lengthsFrom_expected.end (),
                  lengthsFrom_actual1.begin ());
    TEST_EQUALITY( sameLengthsFrom1, true );

    ArrayView<const size_t> lengthsTo_expected = dist_expected.getLengthsTo ();
    ArrayView<const size_t> lengthsTo_actual1 = dist_actual1.getLengthsTo ();
    const bool sameLengthsTo1 =
      lengthsTo_expected.size () == lengthsTo_actual1.size () &&
      std::equal (lengthsTo_expected.begin (),
                  lengthsTo_expected.end (),
                  lengthsTo_actual1.begin ());
    TEST_EQUALITY( sameLengthsTo1, true );

    const bool sameNumReceives2 =
      dist_expected.getNumReceives () == dist_actual2.getNumReceives ();
    TEST_EQUALITY( sameNumReceives2, true );
    const bool sameNumSends2 =
      dist_expected.getNumSends () == dist_actual2.getNumSends ();
    TEST_EQUALITY( sameNumSends2, true );
    const bool sameHasSelfMessage2 =
      dist_expected.hasSelfMessage () == dist_actual2.hasSelfMessage ();
    TEST_EQUALITY( sameHasSelfMessage2, true );
    const bool sameMaxSendLength2 =
      dist_expected.getMaxSendLength () == dist_actual2.getMaxSendLength ();
    TEST_EQUALITY( sameMaxSendLength2, true );
    const bool sameTotalReceiveLength2 =
      dist_expected.getTotalReceiveLength () == dist_actual2.getTotalReceiveLength ();
    TEST_EQUALITY( sameTotalReceiveLength2, true );

    ArrayView<const int> imagesFrom_actual2 = dist_actual2.getImagesFrom ();
    const bool sameImagesFrom2 =
      imagesFrom_expected.size () == imagesFrom_actual2.size () &&
      std::equal (imagesFrom_expected.begin (),
                  imagesFrom_expected.end (),
                  imagesFrom_actual2.begin ());
    TEST_EQUALITY( sameImagesFrom2, true );

    ArrayView<const int> imagesTo_actual2 = dist_actual2.getImagesTo ();
    const bool sameImagesTo2 =
      imagesTo_expected.size () == imagesTo_actual2.size () &&
      std::equal (imagesTo_expected.begin (),
                  imagesTo_expected.end (),
                  imagesTo_actual2.begin ());
    TEST_EQUALITY( sameImagesTo2, true );

    ArrayView<const size_t> lengthsFrom_actual2 = dist_actual2.getLengthsFrom ();
    const bool sameLengthsFrom2 =
      lengthsFrom_expected.size () == lengthsFrom_actual2.size () &&
      std::equal (lengthsFrom_expected.begin (),
                  lengthsFrom_expected.end (),
                  lengthsFrom_actual2.begin ());
    TEST_EQUALITY( sameLengthsFrom2, true );

    ArrayView<const size_t> lengthsTo_actual2 = dist_actual2.getLengthsTo ();
    const bool sameLengthsTo2 =
      lengthsTo_expected.size () == lengthsTo_actual2.size () &&
      std::equal (lengthsTo_expected.begin (),
                  lengthsTo_expected.end (),
                  lengthsTo_actual2.begin ());
    TEST_EQUALITY( sameLengthsTo2, true );
  }
} // namespace (anonymous)

  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP(LOCAL_ORDINAL, GLOBAL_ORDINAL) \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT( ImportUnion, ContigPlusContig, LOCAL_ORDINAL, GLOBAL_ORDINAL )

UNIT_TEST_GROUP(int, int)

UNIT_TEST_GROUP(int, long)

#ifdef HAVE_TEUCHOS_LONG_LONG_INT

// Macros don't like spaces in their arguments.
// typedef long long long_long_type;

// UNIT_TEST_GROUP(int, long_long_type)

#endif // HAVE_TEUCHOS_LONG_LONG_INT



