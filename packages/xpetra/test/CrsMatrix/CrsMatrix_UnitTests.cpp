/*
 * BlockedCrsOperator_UnitTests.cpp
 *
 *  Created on: Aug 22, 2011
 *      Author: wiesner
 */

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Xpetra_ConfigDefs.hpp>

#ifdef HAVE_XPETRA_EPETRAEXT
// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#endif

#include <Xpetra_DefaultPlatform.hpp>
#include <Teuchos_as.hpp>

#include <Xpetra_Map.hpp>
#include <Xpetra_Operator.hpp>
#include <Xpetra_CrsOperator.hpp>
#ifdef HAVE_XPETRA_TPETRA
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_Exceptions.hpp>

#include <MueLu_Utilities.hpp> //TODO: Xpetra tests should not use MueLu

namespace {

  using Teuchos::Array;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::arcp;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Teuchos::Tuple;
  using Teuchos::tuple;
  using std::sort;
  using std::find;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::Comm;

  using Xpetra::DefaultPlatform;
  using Xpetra::Operator;
  using Xpetra::CrsOperator;
#ifdef HAVE_XPETRA_TPETRA
  using Xpetra::TpetraCrsMatrix; //TMP
#endif
  using Xpetra::Map;

  using Xpetra::viewLabel_t;

  bool testMpi = true;
  double errorTolSlack = 1e+1;



  RCP<const Comm<int> > getDefaultComm()
  {
    if (testMpi) {
      return DefaultPlatform::getDefaultPlatform().getComm();
    }
    return rcp(new Teuchos::SerialComm<int>());
  }

  /////////////////////////////////////////////////////

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
                  "test-mpi", "test-serial", &testMpi,
                  "Test MPI (if available) or force test of serial.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
    clp.setOption(
                  "error-tol-slack", &errorTolSlack,
                  "Slack off of machine epsilon used to check test results" );
  }

  //
  // UNIT TESTS
  //


  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, Apply, Scalar, LO, GO, Node )
  {
#ifdef HAVE_XPETRA_EPETRA

    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Operator<Scalar, LO, GO, Node> Operator;
    typedef CrsOperator<Scalar, LO, GO, Node> CrsOperator;
    typedef Xpetra::Map<LO, GO, Node> MapClass;
    typedef Xpetra::MapFactory<LO, GO, Node> MapFactoryClass;

    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();

    Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;

    // generate problem
    LO nEle = 63;
    const RCP<const MapClass> map = MapFactoryClass::Build(lib, nEle, 0, comm);

    RCP<Xpetra::CrsMatrix<Scalar, LO, GO, Node> > matrix =
        Xpetra::CrsMatrixFactory<Scalar,LO,GO,Node>::Build(map, 10);

    LO NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GO> MyGlobalElements = map->getNodeElementList();

    for (LO i = 0; i < NumMyElements; ++i) {
        matrix->insertGlobalValues(MyGlobalElements[i],
                                Teuchos::tuple<GO>(MyGlobalElements[i]),
                                Teuchos::tuple<Scalar>(1.0) );
    }

    matrix->fillComplete();

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(map);

    vec->putScalar(1.0);

    RCP<Xpetra::Vector<Scalar, LO, GO, Node> > vec_sol =
        Xpetra::VectorFactory<Scalar, LO, GO, Node>::Build(matrix->getRangeMap());

    vec_sol->putScalar(0.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, 0.0);

    vec_sol->putScalar(2.0);

    matrix->apply(*vec, *vec_sol, Teuchos::NO_TRANS, 1.0, -0.5);

    TEUCHOS_TEST_COMPARE(vec_sol->norm2(), <, 1e-16, out, success);
#endif
  }




  //
  // INSTANTIATIONS
  //

#define UNIT_TEST_GROUP_ORDINAL( SC, LO, GO, Node )                     \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, Apply, SC, LO, GO, Node )

  typedef Kokkos::DefaultNode::DefaultNodeType DefaultNodeType;

  UNIT_TEST_GROUP_ORDINAL(double, int, int, DefaultNodeType)

}

