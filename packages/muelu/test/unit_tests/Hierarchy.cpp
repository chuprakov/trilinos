#include "Teuchos_UnitTestHarness.hpp"
//#include "Teuchos_ParameterList.hpp"
#include "test_helpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Hierarchy.hpp"
#include "MueLu_PRFactory.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_GenericPRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

TEUCHOS_UNIT_TEST(Hierarchy,Constructor)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Hierarchy> H = rcp(new Hierarchy);

  TEUCHOS_TEST_INEQUALITY(H, Teuchos::null, out, success);

} //Constructor

TEUCHOS_UNIT_TEST(Hierarchy,SetAndGetLevel)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> level = rcp(new Level() );
  H.SetLevel(level);
  RCP<Level> dupLevel = H.GetLevel(0);
  TEUCHOS_TEST_EQUALITY(level, dupLevel, out, success);

}//SetAndGetLevel

TEUCHOS_UNIT_TEST(Hierarchy,NumberOfLevels)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  Hierarchy H;
  RCP<Level> levelOne = rcp(new Level() );
  RCP<Level> levelTwo = rcp(new Level() );
  RCP<Level> levelThree = rcp(new Level() );
  H.SetLevel(levelOne);
  H.SetLevel(levelTwo);
  H.SetLevel(levelThree);
  TEUCHOS_TEST_EQUALITY(H.GetNumberOfLevels(), 3, out, success);

}//NumberOfLevels

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_NoFactoriesGiven)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  out << "Providing no factories to FillHierarchy." << std::endl;
  H.FillHierarchy();

  bool ToF = H.PrintResidualHistory();
  H.PrintResidualHistory(!ToF);
  TEUCHOS_TEST_INEQUALITY(H.PrintResidualHistory(), ToF, out, success);
} //FillHierarchy_NoFactoriesGiven

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_PRFactoryOnly)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);

  out << "Providing just PR factory to FillHierarchy." << std::endl;
  Teuchos::ParameterList status;
  status = H.FillHierarchy(PRFact);
  TEUCHOS_TEST_EQUALITY(status.get("fine nnz",(Cthulhu::global_size_t)-1), 295, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("total nnz",(Cthulhu::global_size_t)-1), 392, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("start level",-1), 0, out, success);
  TEUCHOS_TEST_EQUALITY(status.get("end level",-1), 1, out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(status.get("operator complexity",(SC)-1.0),1.32881,1e-5,out,success);

} //FillHierarchy_PRFactoryOnly

TEUCHOS_UNIT_TEST(Hierarchy,FillHierarchy_BothFactories)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  GenericPRFactory PRFact(PFact);
  RAPFactory    AcFact;

  out << "Providing both factories to FillHierarchy." << std::endl;
  H.FillHierarchy(PRFact,AcFact);
} //FillHierarchy_BothFactories

TEUCHOS_UNIT_TEST(Hierarchy,SetSmoothers)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);
  H.SetSmoothers();

  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPreSmoother()->GetType(),"Ifpack: Gauss-Seidel", out, success);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPostSmoother()->GetType(),"Ifpack: Gauss-Seidel", out, success);


  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Jacobi");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> smooFactory = rcp(new SmootherFactory(smooProto) );
  H.SetSmoothers(*smooFactory);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPreSmoother()->GetType(),"Ifpack: Jacobi", out, success);
  TEUCHOS_TEST_EQUALITY(H.GetLevel(0)->GetPostSmoother()->GetType(),"Ifpack: Jacobi", out, success);

} //SetSmoothers

TEUCHOS_UNIT_TEST(Hierarchy,SetCoarsestSolver)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  SmootherFactory SmooFactory(smoother);

  H.SetCoarsestSolver(SmooFactory);

  RCP<SmootherPrototype>  preSmoo,postSmoo;

  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::PRE);
  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_INEQUALITY(preSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(preSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_EQUALITY(postSmoo, Teuchos::null, out, success);

  H.SetCoarsestSolver(SmooFactory,MueLu::POST);
  preSmoo = levelOne->GetPreSmoother();
  TEUCHOS_TEST_EQUALITY(preSmoo, Teuchos::null, out, success);
  postSmoo = levelOne->GetPostSmoother();
  TEUCHOS_TEST_INEQUALITY(postSmoo, Teuchos::null, out, success);
  TEUCHOS_TEST_EQUALITY(postSmoo->GetType(),"Ifpack: Gauss-Seidel", out, success);

} //SetCoarsestSolver

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_NoArgs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  H.FullPopulate();

} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,FullPopulate_AllArgs)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  RCP<Level> levelOne = rcp(new Level() );
  levelOne->SetLevelID(1);
  RCP<CrsOperator> A = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(99);
  levelOne->SetA(A);

  Hierarchy H;
  H.SetLevel(levelOne);

  RCP<SaPFactory>  PFact = rcp(new SaPFactory());
  RCP<GenericPRFactory> PRFact = rcp(new GenericPRFactory(PFact));
  RCP<RAPFactory>  AcFact = rcp(new RAPFactory());

  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
  RCP<SmootherFactory> SmooFact = rcp( new SmootherFactory(smoother));

  H.FullPopulate(PRFact,AcFact,SmooFact,0,2);

} //FullPopulate

TEUCHOS_UNIT_TEST(Hierarchy,Iterate)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  out << "version: " << MueLu::Version() << std::endl;

  //matrix
  RCP<CrsOperator> Op = MueLu::UnitTest::create_1d_poisson_matrix<SC,LO,GO>(6561);
  RCP<const Map > map = Op->getRowMap();

  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  RCP<Epetra_MultiVector> foo = Utils::MV2NonConstEpetraMV(nullSpace);
  double n;
  foo->Norm1(&n);

  MueLu::Hierarchy<SC,LO,GO,NO,LMO> H;
  H.setDefaultVerbLevel(Teuchos::VERB_HIGH);
  RCP<MueLu::Level<SC,LO,GO,NO,LMO> > Finest = rcp( new MueLu::Level<SC,LO,GO,NO,LMO>() );
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);

  Finest->SetA(Op);
  Finest->Save("Nullspace",nullSpace);
  Finest->Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
                                          //FIXME is implemented

  Finest->Save("NullSpace",nullSpace);
  H.SetLevel(Finest);

  RCP<SaPFactory>         Pfact = rcp( new SaPFactory() );
  RCP<GenericPRFactory>   PRfact = rcp( new GenericPRFactory(Pfact));
  RCP<RAPFactory>         Acfact = rcp( new RAPFactory() );
  Teuchos::ParameterList  ifpackList;
  ifpackList.set("relaxation: type", "Gauss-Seidel");
  ifpackList.set("relaxation: sweeps", (LO) 2);
  ifpackList.set("relaxation: damping factor", (SC) 1.0);
  RCP<SmootherPrototype>  smooProto = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

  RCP<SmootherFactory>    SmooFact = rcp( new SmootherFactory(smooProto) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Teuchos::ParameterList status;
  int maxLevels = 5;
  status = H.FullPopulate(PRfact,Acfact,SmooFact,0,maxLevels);
  out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(out,Teuchos::ParameterList::PrintOptions().indent(2));

  //FIXME we should be able to just call smoother->SetNIts(50) ... but right now an exception gets thrown
  Teuchos::ParameterList amesosList;
  RCP<SmootherPrototype> coarseProto = rcp( new AmesosSmoother("Amesos_Klu",amesosList) );
  SmootherFactory coarseSolveFact(coarseProto);
  H.SetCoarsestSolver(coarseSolveFact,MueLu::PRE);

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
  RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

  RCP<Epetra_MultiVector> epX = Utils::MV2NonConstEpetraMV(X);
  epX->SetSeed(846930886);
  X->randomize();
  //Op->multiply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

  epX->Norm2(&n);
  X->scale(1/n);
  epX->Norm2(&n);
  out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << n << std::endl;

  RHS->putScalar( (SC) 0.0);

  H.PrintResidualHistory(false);
  int iterations=10;
  H.Iterate(*RHS,iterations,*X);

  epX->Norm2(&n);
  out << "||X_" << std::setprecision(2) << iterations << "|| = " << std::setiosflags(ios::fixed) <<
std::setprecision(10) << n << std::endl;

  Teuchos::Array<Teuchos::ScalarTraits<SC>::magnitudeType> norms;
  norms = Utils::ResidualNorm(*Op,*X,*RHS);
  out << "||res_" << std::setprecision(2) << iterations << "|| = " << std::setprecision(15) << norms[0] << std::endl;
  TEUCHOS_TEST_EQUALITY(norms[0]<1e-10, true, out, success);

} //Iterate

}//namespace <anonymous>
