#include "Teuchos_UnitTestHarness.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_Version.hpp"

#include "MueLu_TentativePFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace {

  using Teuchos::rcp;
  using Teuchos::RCP;

  TEUCHOS_UNIT_TEST(TentativePFactory, Constructor)
  {

    out << "version: " << MueLu::Version() << std::endl;

    RCP<TentativePFactory> tentPFact = rcp(new TentativePFactory);
    TEUCHOS_TEST_EQUALITY(tentPFact != Teuchos::null, true, out, success);

  } //Constructor

  TEUCHOS_UNIT_TEST(TentativePFactory, SetGetMethods)
  {
    out << "version: " << MueLu::Version() << std::endl;

    TentativePFactory tentPFact;

    bool flag = tentPFact.TentativeWithQR();
    TEUCHOS_TEST_EQUALITY(flag, false, out, success);
    tentPFact.TentativeWithQR(true);
    flag = tentPFact.TentativeWithQR();
    TEUCHOS_TEST_EQUALITY(flag, true, out, success);
  } //SetGetMethods

  //TODO test BuildP
  //TODO test MakeTentative

  TEUCHOS_UNIT_TEST(TentativePFactory, MakeTentative)
  {
    out << "version: " << MueLu::Version() << std::endl;

    Level fineLevel;
    fineLevel.SetLevelID(1);
    RCP<Operator> A = MueLu::TestHelpers::Factory<SC, LO, GO, NO, LMO>::Build1DPoisson(36);
    fineLevel.SetA(A);

    Level coarseLevel;

    // first iteration calls LAPACK QR
    // second iteration (with only one NS vector) exercises manual orthogonalization
    for (LO NSdim = 2; NSdim >= 1; --NSdim) {
      //RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),1);
      //nullSpace->putScalar( (SC) 1.0);
      // create 3 nullspace vectors
      RCP<MultiVector> nullSpace = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      nullSpace->randomize();
      fineLevel.Save("Nullspace",nullSpace);
      fineLevel.Request("Nullspace"); //FIXME putting this in to avoid error until Merge needs business
      //FIXME is implemented
  
      TentativePFactory::MakeTentative(fineLevel,coarseLevel);

      RCP<Operator> Ptent; 
      coarseLevel.Examine("Ptent",Ptent);
      //Ptent->describe(out,Teuchos::VERB_EXTREME);

      RCP<MultiVector> coarseNullSpace; 
      coarseLevel.Examine("Nullspace",coarseNullSpace);
      //out << *coarseNullSpace << std::endl;

      //check interpolation
      RCP<MultiVector> PtN = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      Ptent->multiply(*coarseNullSpace,*PtN,Teuchos::NO_TRANS,1.0,0.0);

      RCP<MultiVector> diff = MultiVectorFactory::Build(A->getRowMap(),NSdim);
      diff->putScalar(0.0);

      //diff = fineNS + (-1.0)*(P*coarseNS) + 0*diff
      diff->update(1.0,*nullSpace,-1.0,*PtN,0.0);

      //out << *diff << std::endl;

      Teuchos::Array<ST::magnitudeType> norms(NSdim);
      diff->norm2(norms);
      for (LO i=0; i<NSdim; ++i) {
        out << "||diff_" << i << "||_2 = " << norms[i] << std::endl;
        TEUCHOS_TEST_EQUALITY(norms[i]<1e-12, true, out, success);
      }
    } //for (LO NSdim = 1; NSdim <= 2; ++NSdim)

  } //MakeTentative


}//namespace <anonymous>

