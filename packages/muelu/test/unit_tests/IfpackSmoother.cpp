#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_Utilities.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLu::TestHelpers;

  // TODO added by JG: This test should be factorized for every smoothers, to avoid code duplication
  TEUCHOS_UNIT_TEST(IfpackSmoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        GO nx,ny,nz;
        nx = ny = nz = 5;
        GO numGlobalElements = nx*ny*nz;
        RCP<const Map> map = MapFactory::Build(Parameters::getLib(), numGlobalElements, 0, Parameters::getDefaultComm());

        Teuchos::ParameterList  ifpackList;
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );

        RCP<MultiVector> X = MultiVectorFactory::Build(map,1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(map,1);

        //try applying without setting up
        TEST_THROW( smoother->Apply(*X,*RHS) , MueLu::Exceptions::RuntimeError );
        TEST_THROW( smoother->SetNIts(5), MueLu::Exceptions::RuntimeError );
      }
  }

  TEUCHOS_UNIT_TEST(IfpackSmoother, GaussSeidel)
  {
    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {

        //FIXME this will probably fail in parallel b/c it becomes block Jacobi

        out << "version: " << MueLu::Version() << std::endl;
        out << "Tests interface to Ifpack's Gauss-Seidel preconditioner." << std::endl;

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

        Teuchos::ParameterList  ifpackList;
        ifpackList.set("relaxation: type", "Gauss-Seidel");
        ifpackList.set("relaxation: sweeps", (int) 1);
        ifpackList.set("relaxation: damping factor", (double) 1.0);
        ifpackList.set("relaxation: zero starting solution", false);
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("point relaxation stand-alone",ifpackList) );
        Level aLevel;
        aLevel.SetA(Op);

        RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

        smoother->Setup(aLevel);

        X->setSeed(846930886);
        X->randomize();
        Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);

        X->putScalar( (SC) 0.0);

        out << "Applying one GS sweep" << std::endl;
        smoother->Apply(*X,*RHS);
        Teuchos::Array<ST::magnitudeType> norms(1);
        norms = Utils::ResidualNorm(*Op,*X,*RHS);
        TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],6.04555396884098,1e-12,out,success);

        int numIts = 50;
        smoother->SetNIts(numIts);
        TEUCHOS_TEST_EQUALITY(smoother->GetNIts(),50,out,success);
        out << "Applying " << numIts << " GS sweeps" << std::endl;
        X->putScalar( (SC) 0.0);
        smoother->Apply(*X,*RHS);
        norms = Utils::ResidualNorm(*Op,*X,*RHS);
        TEUCHOS_TEST_FLOATING_EQUALITY(norms[0],0.00912675857196253,1e-12,out,success);
      }
  } //GaussSeidel

  TEUCHOS_UNIT_TEST(IfpackSmoother, Chebyshev)
  {

    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {
        out << "version: " << MueLu::Version() << std::endl;

        out << "Tests interface to Ifpack's Chebyshev preconditioner." << std::endl;

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

        Teuchos::ParameterList  ifpackList;
        ifpackList.set("chebyshev: degree", (int) 1);
        ifpackList.set("chebyshev: max eigenvalue", (double) 2.0);
        ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
        ifpackList.set("chebyshev: zero starting solution", false);
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("Chebyshev",ifpackList) );
        Level aLevel;
        aLevel.SetA(Op);

        RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RHS->putScalar( (SC) 0.0);

        smoother->Setup(aLevel);

        X->setSeed(846930886);
        X->randomize();
        Teuchos::Array<ST::magnitudeType> norms(1);
        X->norm2(norms);
        X->scale(1/norms[0]);
        X->norm2(norms);
        out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl; //TODO: why std:ios is on the current namespace ??
        //Op->apply(*X,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);


        int numIts = 1;
        out << "Applying degree " << numIts << " Chebyshev smoother" << std::endl;
        smoother->SetNIts(numIts);
        smoother->Apply(*X,*RHS);
        X->norm2(norms);
        out << "||X_1|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;
        TEUCHOS_TEST_EQUALITY(norms[0]<0.7,true,out,success);  //FIXME should calculate reduction analytically

        numIts = 10;
        smoother->SetNIts(numIts);
        TEUCHOS_TEST_EQUALITY(smoother->GetNIts(),10,out,success);
        out << "Applying degree " << numIts << " Chebyshev smoother" << std::endl;
        X->setSeed(846930886);
        X->randomize();
        X->norm2(norms);
        X->scale(1/norms[0]);
        X->norm2(norms);
        out << "||X_initial|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;
        smoother->Apply(*X,*RHS);
        X->norm2(norms);
        out << "||X_" << std::setprecision(2) << numIts << "|| = " << std::setiosflags(ios::fixed) <<
          std::setprecision(20) << norms[0] << std::endl;
        TEUCHOS_TEST_EQUALITY(norms[0]<0.25,true,out,success);  //FIXME should calculate reduction analytically
      }
  } //Chebyshev

  TEUCHOS_UNIT_TEST(IfpackSmoother, ILU)
  {

    MUELU_TEST_ONLY_FOR(Cthulhu::UseEpetra)
      {
        //FIXME this will probably fail in parallel b/c it becomes block Jacobi

        out << "version: " << MueLu::Version() << std::endl;

        out << "Tests interface to Ifpack's ILU(0) preconditioner." << std::endl;

        RCP<Operator> Op = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(125);

        Teuchos::ParameterList  ifpackList;
        RCP<IfpackSmoother>  smoother = rcp( new IfpackSmoother("ILU",ifpackList) );
        Level aLevel;
        aLevel.SetA(Op);

        RCP<MultiVector> Xtrue = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> X = MultiVectorFactory::Build(Op->getDomainMap(),1);
        RCP<MultiVector> RHS = MultiVectorFactory::Build(Op->getDomainMap(),1);

        smoother->Setup(aLevel);
        Xtrue->setSeed(846930886);
        Xtrue->randomize();
        Teuchos::Array<ST::magnitudeType> norms(1);
        Xtrue->norm2(norms);
        Xtrue->scale(1/norms[0]);
        Op->apply(*Xtrue,*RHS,Teuchos::NO_TRANS,(SC)1.0,(SC)0.0);
        RCP<Epetra_MultiVector> epRHS = Utils::MV2NonConstEpetraMV(RHS);
        RHS->norm2(norms);
        out << "||RHS|| = " << std::setiosflags(ios::fixed) << std::setprecision(10) << norms[0] << std::endl;
        X->putScalar( (SC) 0.0);


        smoother->Apply(*X,*RHS);
        X->norm2(norms);
        out << "||X_final|| = " << std::setiosflags(ios::fixed) << std::setprecision(25) << norms[0] << std::endl;
        norms = Utils::ResidualNorm(*Op,*X,*RHS);
        out << "||residual|| = " << norms[0] << std::endl;
        TEUCHOS_TEST_EQUALITY(norms[0]<1e-10,true,out,success);
        
      }
  } //ILU

}//namespace <anonymous>
