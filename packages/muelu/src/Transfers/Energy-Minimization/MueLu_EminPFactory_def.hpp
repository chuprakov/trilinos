#ifndef MUELU_EMINPFACTORY_DEF_HPP
#define MUELU_EMINPFACTORY_DEF_HPP

#include <Xpetra_CrsGraph_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_EminPFactory_decl.hpp"

#include "MueLu_Utilities.hpp"

#include "MueLu_TentativePFactory.hpp"
#include "MueLu_PatternFactory.hpp"

#include "MueLu_SteepestDescentSolver.hpp"
#include "MueLu_CGSolver.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_Constraint_fwd.hpp"

#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",                 Teuchos::null, "Generating factory for the matrix A used during internal iterations");
    validParamList->set< RCP<const FactoryBase> >("P",                 Teuchos::null, "Generating factory for the initial guess");
    validParamList->set< RCP<const FactoryBase> >("Constraint",        Teuchos::null, "Generating factory for constraints");
    validParamList->set< int >                   ("Niterations",                   3, "Number of iterations of the internal iterative method");
    validParamList->set< int >                   ("Reuse Niterations",             1, "Number of iterations of the internal iterative method");
    validParamList->set< RCP<Matrix> >           ("P0",                Teuchos::null, "Initial guess at P");
    validParamList->set< bool >                  ("Keep P0",                   false, "Keep an initial P0 (for reuse)");
    validParamList->set< RCP<Constraint> >       ("Constraint0",       Teuchos::null, "Initial Constraint");
    validParamList->set< bool >                  ("Keep Constraint0",          false, "Keep an initial Constraint (for reuse)");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");

    static bool isAvailableP0          = false;
    static bool isAvailableConstraint0 = false;

    // Here is a tricky little piece of code
    // We don't want to request (aka call Input) when we reuse and P0 is available
    // However, we cannot run something simple like this:
    //   if (!coarseLevel.IsAvailable("P0", this))
    //     Input(coarseLevel, "P");
    // The reason is that it works fine during the request stage, but fails
    // in the release stage as we _construct_ P0 during Build process. Therefore,
    // we need to understand whether we are in Request or Release mode
    // NOTE: This is a very unique situation, please try not to propagate the
    // mode check any further

    if (coarseLevel.GetRequestMode() == Level::REQUEST) {
      isAvailableP0          = coarseLevel.IsAvailable("P0",          this);
      isAvailableConstraint0 = coarseLevel.IsAvailable("Constraint0", this);
    }

    if (isAvailableP0 == false)
      Input(coarseLevel, "P");

    if (isAvailableConstraint0 == false)
      Input(coarseLevel, "Constraint");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level& coarseLevel) const {
    BuildP(fineLevel, coarseLevel);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EminPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator minimization", coarseLevel);

    const ParameterList & pL = GetParameterList();

    // Get the matrix
    RCP<Matrix> A = Get< RCP<Matrix> >(fineLevel, "A");

    // Get/make initial guess
    RCP<Matrix> P0;
    int         Niterations;
    if (coarseLevel.IsAvailable("P0", this)) {
      // Reuse data
      P0          = coarseLevel.Get<RCP<Matrix> >("P0", this);
      Niterations = pL.get<int>("Reuse Niterations");
      GetOStream(Runtime0, 0) << "Reusing P0" << std::endl;

    } else {
      // Construct data
      P0          = Get< RCP<Matrix> >(coarseLevel, "P");
      Niterations = pL.get<int>("Niterations");
    }
    // NOTE: the main assumption here that P0 satisfies both constraints:
    //   - nonzero pattern
    //   - nullspace preservation

    // Get/make constraint operator
    RCP<Constraint> X;
    if (coarseLevel.IsAvailable("Constraint0", this)) {
      // Reuse data
      X = coarseLevel.Get<RCP<Constraint> >("Constraint0", this);
      GetOStream(Runtime0, 0) << "Reusing Constraint0" << std::endl;

    } else {
      // Construct data
      X = Get< RCP<Constraint> >(coarseLevel, "Constraint");
    }
    GetOStream(Runtime0,0) << "Number of emin iterations = " << Niterations << std::endl;


    RCP<Matrix> P;
    CGSolver EminSolver(Niterations);
    EminSolver.Iterate(*A, *X, *P0, P);

    Set(coarseLevel, "P", P);
    if (pL.get<bool>("Keep P0")) {
      // NOTE: we must do Keep _before_ set as the Needs class only sets if
      //  a) data has been requested (which is not the case here), or
      //  b) data has some keep flag
      coarseLevel.Keep("P0", this);
      Set(coarseLevel, "P0", P);
    }
    if (pL.get<bool>("Keep Constraint0")) {
      // NOTE: we must do Keep _before_ set as the Needs class only sets if
      //  a) data has been requested (which is not the case here), or
      //  b) data has some keep flag
      coarseLevel.Keep("Constraint0", this);
      Set(coarseLevel, "Constraint0", X);
    }

    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    GetOStream(Statistics0,0) << Utils::PrintMatrixInfo(*P, "P", params);
  }

} // namespace MueLu

#endif // MUELU_EMINPFACTORY_DEF_HPP
