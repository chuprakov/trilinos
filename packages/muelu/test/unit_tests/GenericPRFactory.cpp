#include "Teuchos_UnitTestHarness.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

//this macro declares the unit-test-class:
TEUCHOS_UNIT_TEST(GenericPRFactory, Constructor_NoArgs)
{
//we are now in a class method declared by the above macro, and
//that method has these input arguments:
//Teuchos::FancyOStream& out, bool& success

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory());
  TEUCHOS_TEST_EQUALITY(PRFact != Teuchos::null, true, out, success);
} //Constructor_NoArgs

TEUCHOS_UNIT_TEST(GenericPRFactory, Constructor_TwoArgs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<TransPFactory>  RFact = rcp(new TransPFactory());
  RCP<GenericPRFactory>  PRFact = rcp(new GenericPRFactory(PFact,RFact));
  TEUCHOS_TEST_EQUALITY(PRFact != Teuchos::null, true, out, success);
} //Constructor_TwoArgs

TEUCHOS_UNIT_TEST(GenericPRFactory, GetSetMethods)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  GenericPRFactory genericPR = GenericPRFactory();
  genericPR.ReUseAggregates(true);
  TEUCHOS_TEST_EQUALITY( genericPR.ReUseAggregates(), true, out, success);
  genericPR.ReUseAggregates(false);
  TEUCHOS_TEST_EQUALITY( genericPR.ReUseAggregates(), false, out, success);
  genericPR.ReUseGraph(true);
  TEUCHOS_TEST_EQUALITY( genericPR.ReUseGraph(), true, out, success);
  genericPR.ReUseGraph(false);
  TEUCHOS_TEST_EQUALITY( genericPR.ReUseGraph(), false, out, success);
} //GetSetMethods

TEUCHOS_UNIT_TEST(GenericPRFactory, TooCoarse_DoNotBuild)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;
  out << "Test that Build returns early if the coarsest matrix is smaller than specified MaxCoarseSize" << std::endl;

  Level levelOne, levelTwo;
  RCP<Operator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO, GO, Node, LocalMatOps>(240);
  levelOne.SetA(A);

  GenericPRFactory genericPR = GenericPRFactory();
  genericPR.SetMaxCoarseSize(500);
  TEUCHOS_TEST_EQUALITY( genericPR.Build(levelOne,levelTwo), false, out, success);

  out << "Test that Build completes if the coarsest matrix is larger than specified MaxCoarseSize" << std::endl;
  genericPR.SetMaxCoarseSize(239);
  TEUCHOS_TEST_EQUALITY( genericPR.Build(levelOne,levelTwo), true, out, success);
} //TooCoarse_DoNotBuild


}//namespace <anonymous>

