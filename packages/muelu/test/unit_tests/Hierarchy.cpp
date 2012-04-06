#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

TEUCHOS_UNIT_TEST(Hierarchy, Constructor)
{
  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEST_INEQUALITY(H, Teuchos::null);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy, SetAndGetLevel)
{
  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
//   RCP<Level> level = rcp(new Level());
//   H.SetLevel(level);
//   RCP<Level> dupLevel = H.GetLevel(1);

//  TEST_EQUALITY(level, dupLevel);

}//SetAndGetLevel

TEUCHOS_UNIT_TEST(Hierarchy, GetNumLevels)
{

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> levelOne = rcp(new Level() );
  RCP<Level> levelTwo = rcp(new Level() );
  RCP<Level> levelThree = rcp(new Level() );
//   H.SetLevel(levelOne);
//   H.SetLevel(levelTwo);
//   H.SetLevel(levelThree);
//   TEST_EQUALITY(H.GetNumLevels(), 3);

}//GetNumLevels

TEUCHOS_UNIT_TEST(Hierarchy, FillHierarchy_NoFactoriesGiven)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  RCP<Level> & levelOne = H.GetLevel();
  levelOne->Set("A", A);

  out << "Providing no factories to FillHierarchy." << std::endl;
  //  H.FillHierarchy(); TODO

} // FillHierarchy_NoFactoriesGiven


TEUCHOS_UNIT_TEST(Hierarchy, FillHierarchy_BothFactories)
{

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H;
  RCP<Level> levelOne = H.GetLevel();
  levelOne->Set("A", A);

  RCP<SaPFactory>    Pfact = rcp(new SaPFactory());
  RCP<TransPFactory> Rfact = rcp(new TransPFactory(Pfact));
  RAPFactory         AcFact(Pfact, Rfact);

  out << "Providing both factories to FillHierarchy." << std::endl;

  H.FillHierarchy(*Pfact, *Rfact, AcFact);

} //FillHierarchy_BothFactories

TEUCHOS_UNIT_TEST(Hierarchy, SetSmoothers)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H(A);

  RCP<Level> levelTwo = rcp(new Level());
  H.AddLevel(levelTwo);

  // TEST_EQUALITY(H.GetLevel(1)->template Get< RCP<SmootherBase> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel");
  // TEST_EQUALITY(H.GetLevel(1)->template Get< RCP<SmootherBase> >("PostSmoother")->GetType(), "Ifpack: Gauss-Seidel");
  
#ifdef HAVE_MUELU_IFPACK

  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Jacobi");
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto));
  H.SetSmoothers(*smooFactory, 0, 1);
  // JGTODO TEST_EQUALITY(H.GetLevel(1)->Get< RCP<SmootherBase> >("PreSmoother", smooFactory)->GetType(), "Ifpack: Jacobi");
  // JGTODO TEST_EQUALITY(H.GetLevel(1)->Get< RCP<SmootherBase> >("PostSmoother", smooFactory)->GetType(), "Ifpack: Jacobi");
#endif

    }
} //SetSmoothers

TEUCHOS_UNIT_TEST(Hierarchy, SetCoarsestSolver1)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");

  SmootherFactory SmooFactory(smoother);

  Hierarchy H(A);
  RCP<Level> levelOne = H.GetLevel();
  
  //JG FIXME. This is a problem: Users cannot do specific Request(). Even worse if request on two level fact.
  RCP<MueLu::FactoryManagerBase> manager = rcp(new FactoryManager());  // RCP<FactoryManager> manager = rcp(new FactoryManager());
  levelOne->SetFactoryManager(manager);
 
  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  //JG FIXME
  levelOne->SetFactoryManager(Teuchos::null);

  H.SetCoarsestSolver(SmooFactory);

  RCP<SmootherBase>  preSmoo, postSmoo;
  preSmoo = levelOne->Get< RCP<SmootherBase> >("PreSmoother", &SmooFactory);
  TEST_INEQUALITY(preSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(preSmoo->GetType(), "Ifpack: Gauss-Seidel");
  postSmoo = levelOne->Get< RCP<SmootherBase> >("PostSmoother", &SmooFactory);

  TEST_INEQUALITY(postSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(postSmoo->GetType(), "Ifpack: Gauss-Seidel");

#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy, SetCoarsestSolver2)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  SmootherFactory SmooFactory(smoother);

  Hierarchy H(A);
  RCP<Level> levelOne = H.GetLevel();

  //JG FIXME. This is a problem: Users cannot do specific Request(). Even worse if request on two level fact.
  RCP<MueLu::FactoryManagerBase> manager = rcp(new FactoryManager());  // RCP<FactoryManager> manager = rcp(new FactoryManager());
  levelOne->SetFactoryManager(manager);
 
  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  //JG FIXME
  levelOne->SetFactoryManager(Teuchos::null);

  H.SetCoarsestSolver(SmooFactory, MueLu::PRE);

  RCP<SmootherBase> preSmoo;
  preSmoo = levelOne->Get< RCP<SmootherBase> >("PreSmoother", &SmooFactory);

  TEST_INEQUALITY(preSmoo, Teuchos::null);
  //JGTODO  TEST_EQUALITY(preSmoo->GetType(), "Ifpack: Gauss-Seidel");

  TEST_EQUALITY(levelOne->IsAvailable("PostSmoother"), false);

#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy, SetCoarsestSolver3)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  SmootherFactory SmooFactory(smoother);

  Hierarchy H(A);
  RCP<Level> levelOne = H.GetLevel();
  
  //JG FIXME. This is a problem: Users cannot do specific Request(). Even worse if request on two level fact.
  RCP<MueLu::FactoryManagerBase> manager = rcp(new FactoryManager());  // RCP<FactoryManager> manager = rcp(new FactoryManager());
  levelOne->SetFactoryManager(manager);
 
  levelOne->Request("PreSmoother", &SmooFactory);
  levelOne->Request("PostSmoother", &SmooFactory);

  //JG FIXME
  levelOne->SetFactoryManager(Teuchos::null);

  H.SetCoarsestSolver(SmooFactory, MueLu::POST);

  TEST_EQUALITY(levelOne->IsAvailable("PreSmoother"), false);

  RCP<SmootherBase> postSmoo;
  postSmoo = levelOne->Get< RCP<SmootherBase> >("PostSmoother", &SmooFactory);

  TEST_INEQUALITY(postSmoo, Teuchos::null);
  //JG TODO  TEST_EQUALITY(postSmoo->GetType(), "Ifpack: Gauss-Seidel");
#endif
    }
} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy, FullPopulate_NoArgs)
{
  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H(A);

  //TODO  H.FullPopulate();
} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy, FullPopulate_AllArgs)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(99*comm->getSize());

  Hierarchy H(A);

  RCP<SaPFactory>  Pfact = rcp(new SaPFactory());
  RCP<RFactory> Rfact = rcp(new TransPFactory(Pfact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory(Pfact,Rfact));

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));
  H.FullPopulate(*Pfact, *Rfact, *AcFact, *SmooFact, 0, 2);
#endif
    }
} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy, FullPopulate_KeepAggregates)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(399*comm->getSize());

  Hierarchy H(A);
  H.SetMaxCoarseSize(1);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<TentativePFactory>  Pfact = rcp(new TentativePFactory(UCAggFact));
  RCP<RFactory> Rfact = rcp(new TransPFactory(Pfact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory(Pfact,Rfact));

#ifdef HAVE_MUELU_IFPACK
  RCP<SmootherPrototype> smoother = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel");
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));

  H.GetLevel(0)->Keep("Aggregates", UCAggFact.get()); // not working since no FactoryManager available
  H.FullPopulate(*Pfact, *Rfact, *AcFact, *SmooFact, 0, 2);
#endif

  for (LocalOrdinal l=0; l<H.GetNumLevels()-1;l++) {
    TEST_EQUALITY(H.GetLevel(l)->IsAvailable("Aggregates", UCAggFact.get()), true);
  }

    }
} //FullPopulate_KeepAggregates

TEUCHOS_UNIT_TEST(Hierarchy, Iterate)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);

  MueLu::Hierarchy<SC, LO, GO, NO, LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);

  RCP<MueLu::Level> Finest = H.GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->Set("NullSpace", nullSpace);
  Finest->Set("A", Op);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);

  RCP<CoalesceDropFactory> cdFact;
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<RFactory>           Rfact = rcp( new TransPFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory(Pfact,Rfact) );
  

#ifdef HAVE_MUELU_IFPACK
#ifdef HAVE_MUELU_AMESOS
  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel", 2);

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  int maxLevels = 5;
  status = H.FullPopulate(*Pfact, *Rfact, *Acfact, *SmooFact, 0, maxLevels);
  out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(out, Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu", amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact, MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

  X->setSeed(846930886);
  X->randomize();
  //Op->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

  X->norm2(norms);
  X->scale(1/norms[0]);
  X->norm2(norms);
  out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  RHS->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);

  X->norm2(norms);
  out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(std::ios::fixed) <<
    std::setprecision(10) << norms[0] << std::endl;

  norms = Utils::ResidualNorm(*Op, *X, *RHS);
  out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
  TEST_EQUALITY(norms[0]<1e-10, true);

#endif
#endif
    }
} //Iterate

TEUCHOS_UNIT_TEST(Hierarchy, IterateWithImplicitRestriction)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)   //TODO: to be remove in the future
    {

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> Op = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(6561*comm->getSize());  //=8*3^6
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map, 1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);
  nullSpace->norm1(norms);

  MueLu::Hierarchy<SC, LO, GO, NO, LMO> H;
  H.SetImplicitTranspose(true);
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);

  RCP<MueLu::Level> Finest = H.GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A", Op);
  Finest->Set("Nullspace", nullSpace);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  UCAggFact->SetMinNodesPerAggregate(3);
  UCAggFact->SetMaxNeighAlreadySelected(0);
  UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  UCAggFact->SetPhase3AggCreation(0.5);
  RCP<CoalesceDropFactory> cdFact;
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory(TentPFact) );
  RCP<RFactory>           Rfact = rcp( new TransPFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory(Pfact,Rfact) );
  Acfact->SetImplicitTranspose(true);

#ifdef HAVE_MUELU_IFPACK
#ifdef HAVE_MUELU_AMESOS
  RCP<SmootherPrototype> smooProto = TestHelpers::Factory<SC, LO, GO, NO, LMO>::createSmootherPrototype("Gauss-Seidel", 2);

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  int maxLevels = 5;
  status = H.FullPopulate(*Pfact, *Rfact, *Acfact, *SmooFact, 0, maxLevels);
  out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(out, Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu", amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact, MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map, 1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map, 1);

  X->setSeed(846930886);
  X->randomize();
  //Op->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);

  X->norm2(norms);
  X->scale(1/norms[0]);
  X->norm2(norms);
  out << "||X_initial|| = " << std::setiosflags(std::ios::fixed) << std::setprecision(10) << norms[0] << std::endl;

  RHS->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);

  X->norm2(norms);
  out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(std::ios::fixed) <<
    std::setprecision(10) << norms[0] << std::endl;

  norms = Utils::ResidualNorm(*Op, *X, *RHS);
  out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
  TEST_EQUALITY(norms[0]<1e-10, true);

#endif
#endif
    }
} //Iterate

TEUCHOS_UNIT_TEST(Hierarchy, SetupHierarchy1level)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
    {
#ifdef HAVE_MUELU_AMESOS2
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(299*comm->getSize());

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M0; // how to build aggregates and smoother of the first level

  TEST_EQUALITY(H.Setup(0, Teuchos::null,  ptrInArg(M0), Teuchos::null), true);

  RCP<Level> l0 = H.GetLevel(0);

  /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  l0->print(*stdout,Teuchos::VERB_EXTREME);
  l1->print(*stdout,Teuchos::VERB_EXTREME);*/

  TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
  TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

  TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
  TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  //TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

  RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
  RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS->setSeed(846930886);
  RHS->randomize();

  X->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);
#endif
    }
}


TEUCHOS_UNIT_TEST(Hierarchy, SetupHierarchy2level)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
    {
#ifdef HAVE_MUELU_AMESOS
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(299*comm->getSize());

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M0; // how to build aggregates and smoother of the first level

  FactoryManager M1; // first coarse level (Plain aggregation)
  M1.SetFactory("A", rcp(new RAPFactory()));
  M1.SetFactory("P", rcp(new TentativePFactory()));

  FactoryManager M2; // last level (SA)
  M2.SetFactory("A", rcp(new RAPFactory()));
  M2.SetFactory("P", rcp(new SaPFactory()));

  TEST_EQUALITY(H.Setup(0, Teuchos::null,ptrInArg(M0), ptrInArg(M1)), false);
  TEST_EQUALITY(H.Setup(1, ptrInArg(M0), ptrInArg(M1), Teuchos::null), true);

  RCP<Level> l0 = H.GetLevel(0);
  RCP<Level> l1 = H.GetLevel(1);

  /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  l0->print(*stdout,Teuchos::VERB_EXTREME);
  l1->print(*stdout,Teuchos::VERB_EXTREME);*/

  TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
  TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);

  TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
  TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

  TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  // TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

  RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
  RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS->setSeed(846930886);
  RHS->randomize();

  X->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);
#endif
    }
}

TEUCHOS_UNIT_TEST(Hierarchy, SetupHierarchy3level)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
    {
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(299*comm->getSize());

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M0; // how to build aggregates and smoother of the first level

  FactoryManager M1; // first coarse level (Plain aggregation)
  M1.SetFactory("A", rcp(new RAPFactory()));
  RCP<FactoryBase> P = rcp(new TentativePFactory());
  M1.SetFactory("P", P);
  M1.SetFactory("Ptent", P); //FIXME: can it be done automatically in FactoryManager?

  FactoryManager M2; // last level (SA)
  M2.SetFactory("A", rcp(new RAPFactory()));
  M2.SetFactory("P", rcp(new SaPFactory()));

#ifdef HAVE_MUELU_AMESOS2

  TEST_EQUALITY(H.Setup(0, Teuchos::null,  ptrInArg(M0), ptrInArg(M1)), false);
  TEST_EQUALITY(H.Setup(1, ptrInArg(M0), ptrInArg(M1), ptrInArg(M2)), false);
  TEST_EQUALITY(H.Setup(2, ptrInArg(M1), ptrInArg(M2), Teuchos::null ), true);

  RCP<Level> l0 = H.GetLevel(0);
  RCP<Level> l1 = H.GetLevel(1);
  RCP<Level> l2 = H.GetLevel(2);

  /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  l0->print(*stdout,Teuchos::VERB_EXTREME);
  l1->print(*stdout,Teuchos::VERB_EXTREME);
  l2->print(*stdout,Teuchos::VERB_EXTREME);*/

  TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
  TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);

  TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
  TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l0->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

  TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

  TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  // TEST_EQUALITY(l2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

  RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
  RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS->setSeed(846930886);
  RHS->randomize();

  X->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);
#endif
    }
}

TEUCHOS_UNIT_TEST(Hierarchy, SetupHierarchy3levelFacManagers)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
    {
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(299*comm->getSize());

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // setup smoother factory
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 2);
  ifpackList.set("relaxation: damping factor", (SC) 0.9); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  smooProto = Teuchos::rcp( new TrilinosSmoother(ifpackType, ifpackList) );

  RCP<SmootherFactory> preSmooFact;
  RCP<SmootherFactory> postSmooFact;
  preSmooFact = rcp( new SmootherFactory(smooProto) );
  postSmooFact = rcp( new SmootherFactory(smooProto) );
  preSmooFact->SetSmootherPrototypes(smooProto,Teuchos::null );
  postSmooFact->SetSmootherPrototypes(Teuchos::null,smooProto);

  // Multigrid setup phase (using default parameters)
  FactoryManager M0; // how to build aggregates and smoother of the first level
  M0.SetFactory("Smoother", preSmooFact);

  FactoryManager M1; // first coarse level (Plain aggregation)
  M1.SetFactory("A", rcp(new RAPFactory()));
  M1.SetFactory("P", rcp(new TentativePFactory()));
  M1.SetFactory("Smoother", postSmooFact);

  FactoryManager M2; // last level (SA)
  M2.SetFactory("A", rcp(new RAPFactory()));
  M2.SetFactory("P", rcp(new SaPFactory()));

#ifdef HAVE_MUELU_AMESOS
  TEST_EQUALITY(H.Setup(0, Teuchos::null,  ptrInArg(M0), ptrInArg(M1)), false);
  TEST_EQUALITY(H.Setup(1, ptrInArg(M0), ptrInArg(M1), ptrInArg(M2)), false);
  TEST_EQUALITY(H.Setup(2, ptrInArg(M1), ptrInArg(M2), Teuchos::null ), true);

  RCP<Level> l0 = H.GetLevel(0);
  RCP<Level> l1 = H.GetLevel(1);
  RCP<Level> l2 = H.GetLevel(2);

  /*RCP<Teuchos::FancyOStream> stdout = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  l0->print(*stdout,Teuchos::VERB_EXTREME);
  l1->print(*stdout,Teuchos::VERB_EXTREME);
  l2->print(*stdout,Teuchos::VERB_EXTREME);*/

  TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), false);
  TEST_EQUALITY(l2->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false);
  TEST_EQUALITY(l1->IsAvailable("PostSmoother", MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve
  TEST_EQUALITY(l1->IsAvailable("P",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("P",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("R",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("R",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l1->IsAvailable("A",            MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l2->IsAvailable("A",            MueLu::NoFactory::get()), true);

  TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
  TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);

  TEST_EQUALITY(l1->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l1->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final);

  TEST_EQUALITY(l2->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("P",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("R",            MueLu::NoFactory::get()), MueLu::Final);
  TEST_EQUALITY(l2->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  // TEST_EQUALITY(l2->GetKeepFlag("PostSmoother", MueLu::NoFactory::get()), MueLu::Final); // direct solve

  RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
  RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS->setSeed(846930886);
  RHS->randomize();

  X->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);
#endif
    } // test only for Epetra
}

TEUCHOS_UNIT_TEST(Hierarchy, SetupHierarchyTestBreakCondition)
{
  MUELU_TEST_ONLY_FOR(Xpetra::UseEpetra)
    {
  RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  RCP<Operator> A = TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(299*comm->getSize());

  // Multigrid Hierarchy
  Hierarchy H(A);
  H.SetMaxCoarseSize(299*comm->getSize()); // set max coarse size to fit problem size (-> 1 level method)
  H.setVerbLevel(Teuchos::VERB_HIGH);

  // Multigrid setup phase (using default parameters)
  FactoryManager M0; // how to build aggregates and smoother of the first level

  FactoryManager M1; // first coarse level (Plain aggregation)
  M1.SetFactory("A", rcp(new RAPFactory()));
  M1.SetFactory("P", rcp(new TentativePFactory()));

#ifdef HAVE_MUELU_AMESOS
  TEST_EQUALITY(H.Setup(0, Teuchos::null,  ptrInArg(M0), ptrInArg(M1)), true);
  TEST_EQUALITY(H.GetNumLevels(),1);

  RCP<Level> l0 = H.GetLevel(0);
  TEST_EQUALITY(l0->IsAvailable("PreSmoother",  MueLu::NoFactory::get()), true);
  TEST_EQUALITY(l0->IsAvailable("PostSmoother", MueLu::NoFactory::get()), false); // direct solve 
  TEST_EQUALITY(l0->IsAvailable("P",            MueLu::NoFactory::get()), false);
  TEST_EQUALITY(l0->IsAvailable("R",            MueLu::NoFactory::get()), false);
  TEST_EQUALITY(l0->IsAvailable("A",            MueLu::NoFactory::get()), true);

  TEST_EQUALITY(l0->GetKeepFlag("A",            MueLu::NoFactory::get()), MueLu::UserData);
  TEST_EQUALITY(l0->GetKeepFlag("PreSmoother",  MueLu::NoFactory::get()), MueLu::Final);
  // TEST_EQUALITY(l0->GetKeepFlag("PostSmoother",  MueLu::NoFactory::get()), MueLu::Final); //direct solve

  RCP<MultiVector> RHS = MultiVectorFactory::Build(A->getRowMap(), 1);
  RCP<MultiVector> X   = MultiVectorFactory::Build(A->getRowMap(), 1);
  RHS->setSeed(846930886);
  RHS->randomize();

  X->putScalar( (SC) 0.0);

  int iterations=10;
  H.Iterate(*RHS, iterations, *X);
#endif
    } // test only for Epetra
}

}//namespace MueLuTests

//Note from JG:
// For UnitTest,  TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Operator> >("PreSmoother")->GetType(), "Ifpack: Gauss-Seidel");
// should be replaced by
// TEST_EQUALITY(H.GetLevel(1)->Get< RCP<Operator> >("PreSmoother"), preSmoother);
// testing if preSmoother->GetType() == "Ifpack: Gauss-Seidel" should be a unit test of the class IfpackSmoother
