#include <Teuchos_UnitTestHarness.hpp>
#include <Amesos2_config.h>
#include "MueLu_TestHelpers.hpp"
#include "MueLu_TestHelpersSmoothers.hpp"

#include "MueLu_Amesos2Smoother.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

namespace MueLuTests {

  using namespace TestHelpers::Smoothers;

  TEUCHOS_UNIT_TEST(Amesos2Smoother, NotSetup)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
      {
        testApplyNoSetup(Amesos2Smoother(), out, success);
      }
  }

  TEUCHOS_UNIT_TEST(Amesos2Smoother, Apply_Correctness)
  {
    MUELU_TEST_ONLY_FOR(Xpetra::UseTpetra)
      {
#ifdef HAVE_AMESOS2_KLU2
        Amesos2Smoother smoother("Klu");
        testDirectSolver(smoother, out, success);
#endif

#ifdef HAVE_AMESOS2_SUPERLU
        Amesos2Smoother smoother("Superlu");
        testDirectSolver(smoother, out, success);
#endif
      }
  }
  
} // namespace MueLuTests
