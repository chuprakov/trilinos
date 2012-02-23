#include "Teuchos_UnitTestHarness.hpp"

#include "MueLu_config.hpp"

#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Utilities.hpp"
#include "MueLu_PermutedTransferFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Zoltan.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"

#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include "MueLu_GalleryUtils.hpp"

#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {
  
  TEUCHOS_UNIT_TEST(PermutedTransfer, Constructor)
  {
    out << "version: " << MueLu::Version() << std::endl;

    RCP<PermutedTransferFactory> ptFactory = rcp(new PermutedTransferFactory);
    TEST_EQUALITY(ptFactory != Teuchos::null, true);
  } // Constructor test
  
  TEUCHOS_UNIT_TEST(PermutedTransfer, Build1)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel, coarseLevel;
    TestHelpers::Factory<SC, LO, GO, NO, LMO>::createTwoLevelHierarchy(fineLevel, coarseLevel);
    GO nx = 199;
    RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(nx);
    fineLevel.Set("A",A);

    //build coordinates
    Teuchos::ParameterList list;
    list.set("nx",nx);
    RCP<MultiVector> coordVector = MueLu::GalleryUtils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D",A->getRowMap(),list);
    fineLevel.Set("Coordinates",coordVector);


    RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
    RCP<TentativePFactory>    Ptentfact = rcp(new TentativePFactory(UCAggFact));
    RCP<SaPFactory>           Pfact = rcp( new SaPFactory(Ptentfact));
    RCP<RFactory>             Rfact = rcp( new TransPFactory(Pfact) );
    RCP<RAPFactory>           Acfact = rcp( new RAPFactory(Pfact,Rfact) );
    RCP<RFactory>             Rtentfact = rcp( new TransPFactory(Ptentfact) );

    RCP<MultiVectorTransferFactory> mvTransFact = rcp(new MultiVectorTransferFactory("Coordinates","R",Rtentfact));
    Acfact->AddTransferFactory(mvTransFact);
    RCP<ZoltanInterface>      zoltan = rcp(new ZoltanInterface(Acfact,mvTransFact));
    RCP<RepartitionFactory> RepartitionFact = rcp(new RepartitionFactory(zoltan,Acfact));

    coarseLevel.Request("A",Acfact.get());  // kick off the DeclareInputs
    coarseLevel.Request("Permutation",RepartitionFact.get());  // request permutation matrix
    //coarseLevel.Request("P",Pfact.get());
    coarseLevel.Request("R",Rtentfact.get());
    coarseLevel.Request("Coordinates",mvTransFact.get());

    RCP<PermutedTransferFactory> ptFactory = rcp( new PermutedTransferFactory(RepartitionFact, Acfact, Pfact, MueLu::INTERPOLATION) );
    coarseLevel.Request("P",ptFactory.get());
    ptFactory->Build(fineLevel,coarseLevel);

    ptFactory = rcp( new PermutedTransferFactory(RepartitionFact, Acfact, Rfact, MueLu::RESTRICTION) );
    coarseLevel.Request("R",ptFactory.get());
    ptFactory->Build(fineLevel,coarseLevel);



  } // Constructor test

} // namespace MueLuTests

