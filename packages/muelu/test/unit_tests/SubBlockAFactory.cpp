#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>

#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  /////////////////////////
  // helper function

  Teuchos::RCP<CrsOperator> GenerateProblemMatrix(const Teuchos::RCP<const Map> map, Scalar a = 2.0, Scalar b = -1.0, Scalar c = -1.0) {


    Teuchos::RCP<CrsOperator> mtx = MueLu::Gallery::MatrixTraits<Map,CrsOperator>::Build(map, 3);

    LocalOrdinal NumMyElements = map->getNodeNumElements();
    Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
    GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
    GlobalOrdinal nIndexBase = map->getIndexBase();

    GlobalOrdinal NumEntries;
    LocalOrdinal nnz=2;
    std::vector<Scalar> Values(nnz);
    std::vector<GlobalOrdinal> Indices(nnz);

    for (LocalOrdinal i = 0; i < NumMyElements; ++i)
    {
      if (MyGlobalElements[i] == nIndexBase)
      {
        // off-diagonal for first row
        Indices[0] = nIndexBase;
        NumEntries = 1;
        Values[0] = c;
      }
      else if (MyGlobalElements[i] == nIndexBase + NumGlobalElements - 1)
      {
        // off-diagonal for last row
        Indices[0] = nIndexBase + NumGlobalElements - 2;
        NumEntries = 1;
        Values[0] = b;
      }
      else
      {
        // off-diagonal for internal row
        Indices[0] = MyGlobalElements[i] - 1;
        Values[1] = b;
        Indices[1] = MyGlobalElements[i] + 1;
        Values[0] = c;
        NumEntries = 2;
      }

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

      // Put in the diagonal entry
      mtx->insertGlobalValues(MyGlobalElements[i],
          Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
          Teuchos::tuple<Scalar>(a) );

    } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)

    mtx->fillComplete(map,map);
    
    return mtx;
  }

  TEUCHOS_UNIT_TEST(SubBlockAFactory, Constructor)
  {
    // test for accessing subblocks from a blocked CRS Operator using SubBlockAFactory
    out << "version: " << MueLu::Version() << std::endl;

    RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

    /**********************************************************************************/
    /* CREATE INITIAL MATRIX                                                          */
    /**********************************************************************************/
    RCP<const Map> bigMap;
    RCP<const Map> map1;
    RCP<const Map> map2;
    GO numElements = 500;
    GO numElements1 = 400;
    GO numElements2 = 100;

    Xpetra::UnderlyingLib lib = MueLuTests::TestHelpers::Parameters::getLib();
    
    map1   = MapFactory::Build(lib, numElements1, 0, comm);
    map2   = MapFactory::Build(lib, numElements2, numElements1, comm);

    std::vector<GlobalOrdinal> localGids; // vector with all local GIDs on cur proc
    Teuchos::ArrayView< const GlobalOrdinal > map1eleList = map1->getNodeElementList(); // append all local gids from map1 and map2
    localGids.insert(localGids.end(), map1eleList.begin(), map1eleList.end());
    Teuchos::ArrayView< const GlobalOrdinal > map2eleList = map2->getNodeElementList();
    localGids.insert(localGids.end(), map2eleList.begin(), map2eleList.end());
    Teuchos::ArrayView<GlobalOrdinal> eleList(&localGids[0],localGids.size());
    bigMap = MapFactory::Build(lib, numElements, eleList, 0, comm); // create full big map (concatenation of map1 and map2)

    std::vector<Teuchos::RCP<const Map> > maps;
    maps.push_back(map1); maps.push_back(map2);

    Teuchos::RCP<const Xpetra::MapExtractor<Scalar, LO, GO, Node> > mapExtractor = Xpetra::MapExtractorFactory<Scalar,LO,GO,Node>::Build(bigMap, maps);

    RCP<CrsOperator> Op11 = GenerateProblemMatrix(map1,2,-1,-1);
    RCP<CrsOperator> Op22 = GenerateProblemMatrix(map2,3,-2,-1);
    
    // build blocked operator
    Teuchos::RCP<Xpetra::BlockedCrsOperator<Scalar,LO,GO,Node,LocalMatOps> > bOp = Teuchos::rcp(new Xpetra::BlockedCrsOperator<Scalar,LO,GO,Node,LocalMatOps>(mapExtractor,mapExtractor,10));

    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat11 = Op11->getCrsMatrix();
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LO,GO,Node,LocalMatOps> > crsMat22 = Op22->getCrsMatrix();
    bOp->setMatrix(0,0,crsMat11);
    bOp->setMatrix(1,1,crsMat22);
    bOp->fillComplete();
    TEST_EQUALITY(bOp!=Teuchos::null, true);

    // build hierarchy
    RCP<Level> levelOne = rcp(new Level());
    levelOne->Set("A", Teuchos::rcp_dynamic_cast<Operator>(bOp)); // set blocked operator

    // define sub block factories for blocked operator "A"
    RCP<SubBlockAFactory> A11Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 0, 0));
    RCP<SubBlockAFactory> A22Fact = Teuchos::rcp(new SubBlockAFactory(MueLu::NoFactory::getRCP(), 1, 1));

    // request subblocks of A
    levelOne->Request("A", A11Fact.get(), MueLu::NoFactory::get());
    levelOne->Request("A", A22Fact.get(), MueLu::NoFactory::get());
    TEST_EQUALITY(levelOne->IsRequested("A", A11Fact.get()),true);
    TEST_EQUALITY(levelOne->IsRequested("A", A22Fact.get()),true);

    RCP<Operator> A11 = levelOne->Get<RCP<Operator> >("A",A11Fact.get());
    RCP<Operator> A22 = levelOne->Get<RCP<Operator> >("A",A22Fact.get());
    TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()),true);
    TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()),true);

    levelOne->Release("A", A11Fact.get());
    levelOne->Release("A", A22Fact.get());

    TEST_EQUALITY(levelOne->IsAvailable("A", A11Fact.get()),false);
    TEST_EQUALITY(levelOne->IsAvailable("A", A22Fact.get()),false);
    TEST_EQUALITY(A11->getRowMap()->isSameAs(*(Op11->getRowMap())), true);
    TEST_EQUALITY(A11->getColMap()->isSameAs(*(Op11->getColMap())), true);
    TEST_EQUALITY(A11->getRangeMap()->isSameAs(*(Op11->getRangeMap())), true);
    TEST_EQUALITY(A11->getDomainMap()->isSameAs(*(Op11->getDomainMap())), true);
    TEST_EQUALITY(A11->getNodeNumEntries(),Op11->getNodeNumEntries());
    TEST_EQUALITY(A22->getRowMap()->isSameAs(*(Op22->getRowMap())), true);
    TEST_EQUALITY(A22->getColMap()->isSameAs(*(Op22->getColMap())), true);
    TEST_EQUALITY(A22->getRangeMap()->isSameAs(*(Op22->getRangeMap())), true);
    TEST_EQUALITY(A22->getDomainMap()->isSameAs(*(Op22->getDomainMap())), true);
    TEST_EQUALITY(A22->getNodeNumEntries(),Op22->getNodeNumEntries());

    //levelOne->print(out,Teuchos::VERB_EXTREME);
  } //Constructor
}


