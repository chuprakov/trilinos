// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP
#define MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP

#include <Teuchos_XMLParameterListHelpers.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_EasyParameterListInterpreter_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_FactoryManager.hpp"

#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoarseMapFactory.hpp"
#include "MueLu_CoordinatesTransferFactory.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_FilteredAFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_RebalanceAcFactory.hpp"
#include "MueLu_RebalanceTransferFactory.hpp"
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_Zoltan2Interface.hpp"

#ifdef HAVE_MUELU_EXPERIMENTAL
#include "MueLu_ConstraintFactory.hpp"
#include "MueLu_PatternFactory.hpp"
#include "MueLu_EminPFactory.hpp"
#endif

namespace MueLu {

  //! Helper functions to compare two paramter lists
  static bool areSame(const ParameterList& list1, const ParameterList& list2);


  // This macro is tricky. The use case is when we do not have a level specific parameter, so we
  // need to take the default value from the general list for all levels, or set it to the default.
#define MUELU_READ_2LIST_PARAM(paramList, defaultList, paramStr, varType, defaultValue, varName) \
  varType varName; \
  if      (paramList.isParameter(paramStr))   varName = paramList.get<varType>(paramStr); \
  else if (defaultList.isParameter(paramStr)) varName = defaultList.get<varType>(paramStr); \
  else                                        varName = paramList.get<varType>(paramStr, defaultValue);

  // This macro check whether the variable is in the list.
  // If it is, it copies its value to the second list, possibly with a new name
  // Similar to the above macro, we all try to take a value from the default list
  // NOTE: this essentially converts UserAPI parameter names into MueLu internal ones
#define MUELU_TEST_AND_SET_PARAM(listWrite, varNameWrite, paramList, defaultList, varNameRead, T) \
  if      (paramList.isParameter(varNameRead))   listWrite.set(varNameWrite, paramList.get<T>(varNameRead)); \
  else if (defaultList.isParameter(varNameRead)) listWrite.set(varNameWrite, defaultList.get<T>(varNameRead));

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EasyParameterListInterpreter(Teuchos::ParameterList& paramList) {
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::EasyParameterListInterpreter(const std::string& xmlFileName, const Teuchos::Comm<int>& comm) {
    Teuchos::ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<Teuchos::ParameterList>(&paramList), comm);
    SetParameterList(paramList);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameterList(const Teuchos::ParameterList& paramList_) {
    // Create a non const copy of the parameter list
    // Working with a modifiable list is much much easier than with original one
    ParameterList paramList = paramList_;

    // Translate cycle type parameter
    Cycle_ = Hierarchy::GetDefaultCycle();
    if (paramList.isParameter("cycle type")) {
      std::map<std::string,CycleType> cycleMap;
      cycleMap["V"] = VCYCLE;
      cycleMap["W"] = WCYCLE;

      std::string cycleType = paramList.get<std::string>("cycle type");
      TEUCHOS_TEST_FOR_EXCEPTION(cycleMap.count(cycleType) == 0, Exceptions::RuntimeError, "Invalid cycle type: \"" << cycleType << "\"");
      Cycle_ = cycleMap[cycleType];
    }

    this->maxCoarseSize_       = paramList.get<int> ("coarse: max size",    Hierarchy::GetDefaultMaxCoarseSize());
    this->numDesiredLevel_     = paramList.get<int> ("max levels",          Hierarchy::GetDefaultMaxLevels());
    this->graphOutputLevel_    = paramList.get<int> ("debug: graph level", -1);
    this->blockSize_           = paramList.get<int> ("number of equations", 1);

    // Save level data
    if (paramList.isSublist("print")) {
      ParameterList printList = paramList.sublist("print");

      if (printList.isParameter("A"))
        this->matricesToPrint_     = Teuchos::getArrayFromStringParameter<int>(printList, "A");
      if (printList.isParameter("P"))
        this->prolongatorsToPrint_ = Teuchos::getArrayFromStringParameter<int>(printList, "P");
      if (printList.isParameter("R"))
        this->restrictorsToPrint_  = Teuchos::getArrayFromStringParameter<int>(printList, "R");
    }

    // Translate verbosity parameter
    this->verbosity_ = static_cast<MsgType>(Hierarchy::GetDefaultVerbLevel());      // cast int to enum
    if (paramList.isParameter("verbosity")) {
      std::map<std::string,MsgType> verbMap;
      verbMap["none"]    = None;
      verbMap["low"]     = Low;
      verbMap["medium"]  = Medium;
      verbMap["high"]    = High;
      verbMap["extreme"] = Extreme;
      verbMap["test"]    = Test;

      std::string verbosityLevel = paramList.get<std::string>("verbosity");
      TEUCHOS_TEST_FOR_EXCEPTION(verbMap.count(verbosityLevel) == 0, Exceptions::RuntimeError, "Invalid verbosity level: \"" << verbosityLevel << "\"");
      this->verbosity_ = verbMap[verbosityLevel];
    }

    // Detect if we need to transfer coordinates to coarse levels. We do that iff
    //  - we use "laplacian" dropping on some level, or
    //  - we use repartitioning on some level
    // This is not ideal, as we may have "repartition: enable" turned on by default
    // and not present in the list, but it is better than nothing.
    useCoordinates_ = false;
    if ((paramList.isParameter("repartition: enable")      && paramList.get<bool>("repartition: enable")             == true) ||
        (paramList.isParameter("aggregation: drop scheme") && paramList.get<std::string>("aggregation: drop scheme") == "laplacian")) {
      useCoordinates_ = true;

    } else {
      for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
        std::string levelStr = "level" + toString(levelID);

        if (paramList.isSublist(levelStr)) {
          const ParameterList& levelList = paramList.sublist(levelStr);

          if ((levelList.isParameter("repartition: enable")      && levelList.get<bool>("repartition: enable")             == true) ||
              (levelList.isParameter("aggregation: drop scheme") && levelList.get<std::string>("aggregation: drop scheme") == "laplacian")) {
            useCoordinates_ = true;
            break;
          }
        }
      }
    }

    // Detect if we do implicit P and R rebalance
    if (paramList.isParameter("repartition: enable") && paramList.get<bool>("repartition: enable") == true)
      this->implicitPRrebalance_ = !paramList.get<bool>("repartition: rebalance P and R", Hierarchy::GetDefaultPRrebalance());

    // Create default manager
    RCP<FactoryManager> defaultManager = rcp(new FactoryManager());
    defaultManager->SetVerbLevel(this->verbosity_);
    UpdateFactoryManager(paramList, ParameterList(), *defaultManager);
    defaultManager->Print();

    for (int levelID = 0; levelID < this->numDesiredLevel_; levelID++) {
      RCP<FactoryManager> levelManager;

      if (paramList.isSublist("level " + toString(levelID))) {
        // Some level specific parameters, update default manager
        bool mustAlreadyExist = true;
        ParameterList& levelList = paramList.sublist("level " + toString(levelID), mustAlreadyExist);

        levelManager = rcp(new FactoryManager(*defaultManager));
        levelManager->SetVerbLevel(defaultManager->GetVerbLevel());

        UpdateFactoryManager(levelList, paramList, *levelManager);

      } else {
        // No level specific parameter, use default manager
        levelManager = defaultManager;
      }

      this->AddFactoryManager(levelID, 1, levelManager);
    }
    // FIXME: parameters passed to packages, like Ifpack2, are not touched by us, resulting in "[unused]" flag
    // being displayed. On the other hand, we don't want to simply iterate through them touching. I don't know
    // what a good solution looks like
    this->GetOStream(static_cast<MsgType>(Runtime1 | Test), 0) << paramList << std::endl;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupMatrix(Matrix& A) const {
    if (A.GetFixedBlockSize() != blockSize_)
      this->GetOStream(Warnings0) << "Warning: setting matrix block size to " << blockSize_ << " (value of \"number of equations\" parameter in the list) "
          << "instead of " << A.GetFixedBlockSize() << " (provided matrix)." << std::endl;
    A.SetFixedBlockSize(blockSize_);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetupHierarchy(Hierarchy& H) const {
    HierarchyManager::SetupHierarchy(H);
    H.SetCycle(Cycle_);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void EasyParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::UpdateFactoryManager(Teuchos::ParameterList& paramList,
        const Teuchos::ParameterList& defaultList, FactoryManager& manager) {
    // NOTE: Factory::SetParameterList must be called prior to Factory::SetFactory, as
    // SetParameterList sets default values for non mentioned parameters, including factories

    // === Smoothing ===
    bool isCustomSmoother =
        paramList.isParameter("smoother: pre or post") ||
        paramList.isParameter("smoother: type")   || paramList.isParameter("smoother: pre type")   || paramList.isParameter("smoother: post type") ||
        paramList.isSublist("smoother: params")   || paramList.isSublist("smoother: pre params")   || paramList.isSublist("smoother: post params") ||
        paramList.isParameter("smoother: sweeps") || paramList.isParameter("smoother: pre sweeps") || paramList.isParameter("smoother: post sweeps");
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: pre or post", std::string, "both", PreOrPost);
    if (PreOrPost == "none") {
      manager.SetFactory("Smoother", Teuchos::null);

    } else if (isCustomSmoother) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("smoother: type") && paramList.isParameter("smoother: pre type"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: type\" and \"smoother: pre type\"");
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("smoother: type") && paramList.isParameter("smoother: post type"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: type\" and \"smoother: post type\"");
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("smoother: sweeps") && paramList.isParameter("smoother: pre sweeps"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: sweeps\" and \"smoother: pre sweeps\"");
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isParameter("smoother: sweeps") && paramList.isParameter("smoother: post sweeps"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: sweeps\" and \"smoother: post sweeps\"");
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: params") && paramList.isSublist("smoother: pre params"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: params\" and \"smoother: pre params\"");
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.isSublist("smoother: params") && paramList.isSublist("smoother: pre params"),
                                 Exceptions::InvalidArgument, "You cannot specify both \"smoother: params\" and \"smoother: pre params\"");
      TEUCHOS_TEST_FOR_EXCEPTION(PreOrPost == "both" && (paramList.isParameter("smoother: pre type") != paramList.isParameter("smoother: post type")),
                                 Exceptions::InvalidArgument, "You must specify both \"smoother: pre type\" and \"smoother: post type\"");

      // Default values
      const int overlap = 0;
      ParameterList defaultSmootherParams;
      defaultSmootherParams.set("relaxation: type",           "Symmetric Gauss-Seidel");
      defaultSmootherParams.set("relaxation: sweeps",         Teuchos::OrdinalTraits<LO>::one());
      defaultSmootherParams.set("relaxation: damping factor", Teuchos::ScalarTraits<Scalar>::one());

      RCP<SmootherPrototype> preSmoother = Teuchos::null, postSmoother = Teuchos::null;
      std::string            preSmootherType,             postSmootherType;
      ParameterList          preSmootherParams,           postSmootherParams;
      if (PreOrPost == "pre" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: pre type")) {
          preSmootherType = paramList.get<std::string>("smoother: pre type");
        } else {
          MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: type", std::string, "RELAXATION", preSmootherTypeTmp);
          preSmootherType = preSmootherTypeTmp;
        }

        if (paramList.isSublist("smoother: pre params"))
          preSmootherParams = paramList.sublist("smoother: pre params");
        else if (paramList.isSublist("smoother: params"))
          preSmootherParams = paramList.sublist("smoother: params");
        else if (defaultList.isSublist("smoother: params"))
          preSmootherParams = defaultList.sublist("smoother: params");
        else if (preSmootherType == "RELAXATION")
          preSmootherParams = defaultSmootherParams;

        preSmoother = rcp(new TrilinosSmoother(preSmootherType, preSmootherParams, overlap));
      }

      if (PreOrPost == "post" || PreOrPost == "both") {
        if (paramList.isParameter("smoother: post type"))
          postSmootherType = paramList.get<std::string>("smoother: post type");
        else {
          MUELU_READ_2LIST_PARAM(paramList, defaultList, "smoother: type", std::string, "RELAXATION", postSmootherTypeTmp);
          postSmootherType = postSmootherTypeTmp;
        }

        if (paramList.isSublist("smoother: post params"))
          postSmootherParams = paramList.sublist("smoother: post params");
        else if (paramList.isSublist("smoother: params"))
          postSmootherParams = paramList.sublist("smoother: params");
        else if (defaultList.isSublist("smoother: params"))
          postSmootherParams = defaultList.sublist("smoother: params");
        else if (postSmootherType == "RELAXATION")
          postSmootherParams = defaultSmootherParams;

        if (postSmootherType == preSmootherType && areSame(preSmootherParams, postSmootherParams))
          postSmoother = preSmoother;
        else
          postSmoother = rcp(new TrilinosSmoother(postSmootherType, postSmootherParams, overlap));
      }

      manager.SetFactory("Smoother", rcp(new SmootherFactory(preSmoother, postSmoother)));
    }

    // === Coarse solver ===
    bool isCustomCoarseSolver =
        paramList.isParameter("coarse: type")   ||
        paramList.isParameter("coarse: params");
    if (paramList.isParameter("coarse: type") && paramList.get<std::string>("coarse: type") == "none") {
      manager.SetFactory("CoarseSolver", Teuchos::null);

    } else if (isCustomCoarseSolver) {
      // FIXME: get default values from the factory
      // NOTE: none of the smoothers at the moment use parameter validation framework, so we
      // cannot get the default values from it.
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "coarse: type", std::string, "", coarseType);

      ParameterList coarseParams;
      if (paramList.isSublist("coarse: params"))
        coarseParams = paramList.sublist("coarse: params");
      else if (defaultList.isSublist("coarse: params"))
        coarseParams = defaultList.sublist("coarse: params");

      RCP<SmootherPrototype> coarseSmoother;
      // TODO: this is not a proper place to check. If we consider direct solver to be a special
      // case of smoother, we would like to unify Amesos and Ifpack2 smoothers in src/Smoothers, and
      // have a single factory responsible for those. Then, this check would belong there.
      if (coarseType == "RELAXATION" || coarseType == "CHEBYSHEV" ||
          coarseType == "ILUT" || coarseType == "ILU" || coarseType == "RILUK")
        coarseSmoother = rcp(new TrilinosSmoother(coarseType, coarseParams));
      else
        coarseSmoother = rcp(new DirectSolver(coarseType, coarseParams));

      manager.SetFactory("CoarseSolver", rcp(new SmootherFactory(coarseSmoother)));
    }

    // === Aggregation ===
    // Aggregation graph
    RCP<CoalesceDropFactory> dropFactory = rcp(new CoalesceDropFactory());
    ParameterList dropParams = *(dropFactory->GetValidParameterList());
    dropParams.set                      ("lightweight wrap", true);
    MUELU_TEST_AND_SET_PARAM(dropParams, "algorithm",                     paramList, defaultList, "aggregation: drop scheme",         std::string);
    // Rename classical to original
    if (dropParams.isParameter("algorithm") && dropParams.get<std::string>("algorithm") == "classical")
      dropParams.set("algorithm", "original");
    MUELU_TEST_AND_SET_PARAM(dropParams, "aggregation threshold",         paramList, defaultList, "aggregation: drop tol",            double);
    MUELU_TEST_AND_SET_PARAM(dropParams, "Dirichlet detection threshold", paramList, defaultList, "aggregation: Dirichlet threshold", double);

    dropFactory->SetParameterList(dropParams);
    manager.SetFactory("Graph", dropFactory);

    // Aggregation sheme
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "aggregation: type", std::string, "uncoupled", aggType);
    RCP<Factory> aggFactory;
    if      (aggType == "uncoupled") aggFactory = rcp(new UncoupledAggregationFactory());
    else if (aggType == "coupled")   aggFactory = rcp(new CoupledAggregationFactory());
    aggFactory->SetFactory("Graph",       manager.GetFactory("Graph"));
    aggFactory->SetFactory("DofsPerNode", manager.GetFactory("Graph"));
    manager.SetFactory("Aggregates", aggFactory);

    // Coarse map
    RCP<CoarseMapFactory> coarseMap = rcp(new CoarseMapFactory());
    coarseMap->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    manager.SetFactory("CoarseMap", coarseMap);

    // Tentative P
    RCP<TentativePFactory> Ptent = rcp(new TentativePFactory());
    Ptent->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
    Ptent->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
    manager.SetFactory("Ptent",     Ptent);

    // Nullspace
    RCP<NullspaceFactory> nullSpace = rcp(new NullspaceFactory());
    nullSpace->SetFactory("Nullspace", manager.GetFactory("Ptent"));
    manager.SetFactory("Nullspace", nullSpace);

    // === Prolongation ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "multigrid algorithm", std::string, "sa", multigridAlgo);
    if (multigridAlgo == "sa") {
      // Smoothed aggregation
      RCP<SaPFactory> P = rcp(new SaPFactory());
      ParameterList Pparams = *(P->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(Pparams, "Damping factor", paramList, defaultList, "sa: damping factor", double);
      P->SetParameterList(Pparams);

      if (paramList.isParameter("sa: use filtered matrix") && paramList.get<bool>("sa: use filtered matrix")) {
        // Filtering
        RCP<FilteredAFactory> filterFactory = rcp(new FilteredAFactory());
        ParameterList fParams = *(filterFactory->GetValidParameterList());
        MUELU_TEST_AND_SET_PARAM(fParams, "lumping", paramList, defaultList, "filtered matrix: use lumping", bool);
        filterFactory->SetParameterList(fParams);
        filterFactory->SetFactory("Graph", manager.GetFactory("Graph"));
        P->SetFactory("A", filterFactory);
      }

      P->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("P", P);

    } else if (multigridAlgo == "emin") {
#ifdef HAVE_MUELU_EXPERIMENTAL
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "emin: pattern", std::string, "AkPtent", patternType);
      TEUCHOS_TEST_FOR_EXCEPTION(patternType != "AkPtent", Exceptions::InvalidArgument, "Invalid pattern name: \"" << patternType << "\". Valid options: \"AkPtent\"");

      // Pattern
      RCP<PatternFactory> patternFactory = rcp(new PatternFactory());
      ParameterList patternParams = *(patternFactory->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(patternParams, "k", paramList, defaultList, "emin: pattern order", int);
      patternFactory->SetParameterList(patternParams);
      patternFactory->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("Ppattern", patternFactory);

      // Constraint
      RCP<ConstraintFactory> constraintFactory = rcp(new ConstraintFactory());
      constraintFactory->SetFactory("Ppattern",        manager.GetFactory("Ppattern"));
      constraintFactory->SetFactory("CoarseNullspace", manager.GetFactory("Ptent"));
      manager.SetFactory("Constraint", constraintFactory);

      // Energy minimization
      RCP<EminPFactory> P = rcp(new EminPFactory());
      ParameterList Pparams = *(P->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(Pparams, "Niterations", paramList, defaultList, "emin: num iterations", int);
      P->SetParameterList(Pparams);
      P->SetFactory("P",          manager.GetFactory("Ptent"));
      P->SetFactory("Constraint", manager.GetFactory("Constraint"));
      manager.SetFactory("P", P);
#else
      throw Exceptions::RuntimeError("Please enable Experimental options in MueLu to use \"emin\"");
#endif

    } else if (multigridAlgo == "pg") {
      // Petrov-Galerkin
      RCP<PgPFactory> P = rcp(new PgPFactory());
      P->SetFactory("P", manager.GetFactory("Ptent"));
      manager.SetFactory("P", P);
    }

    // === Restriction ===
    RCP<TransPFactory> R = rcp(new TransPFactory());
    R->SetFactory("P", manager.GetFactory("P"));
    manager.SetFactory("R", R);

    // === RAP ===
    RCP<RAPFactory> RAP = rcp(new RAPFactory());
    RAP->SetFactory("P", manager.GetFactory("P"));
    RAP->SetFactory("R", manager.GetFactory("R"));
    manager.SetFactory("A", RAP);

    // === Coordinates ===
    if (useCoordinates_) {
      RCP<CoordinatesTransferFactory> coords = rcp(new CoordinatesTransferFactory());
      coords->SetFactory("Aggregates", manager.GetFactory("Aggregates"));
      coords->SetFactory("CoarseMap",  manager.GetFactory("CoarseMap"));
      manager.SetFactory("Coordinates", coords);

      RAP->AddTransferFactory(manager.GetFactory("Coordinates"));
    }

    // === Repartitioning ===
    MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: enable", bool, false, enableRepart);
    if (enableRepart) {
#ifdef HAVE_MPI
      MUELU_READ_2LIST_PARAM(paramList, defaultList, "repartition: partitioner", std::string, "zoltan", partName);
      TEUCHOS_TEST_FOR_EXCEPTION(partName != "zoltan" && partName != "zoltan2", Exceptions::InvalidArgument,
                                 "Invalid partitioner name: \"" << partName << "\". Valid options: \"zoltan\", \"zoltan2\"");
      // Partitioner
      RCP<Factory> partitioner;
      if (partName == "zoltan") {
#ifdef HAVE_MUELU_ZOLTAN
        partitioner = rcp(new ZoltanInterface());
        // NOTE: ZoltanInteface ("zoltan") does not support external parameters through ParameterList
#else
        throw Exceptions::RuntimeError("Zoltan interface is not available");
#endif
      } else if (partName == "zoltan2") {
#ifdef HAVE_MUELU_ZOLTAN2
        partitioner = rcp(new Zoltan2Interface());
        ParameterList partParams = *(partitioner->GetValidParameterList());
        RCP<const ParameterList> partpartParams = rcp(new ParameterList(paramList.sublist("repartition: params", false)));
        partParams.set("ParameterList", partpartParams);
        partitioner->SetParameterList(partParams);
#else
        throw Exceptions::RuntimeError("Zoltan2 interface is not available");
#endif
      }
      partitioner->SetFactory("A",           manager.GetFactory("A"));
      partitioner->SetFactory("Coordinates", manager.GetFactory("Coordinates"));
      manager.SetFactory("Partition", partitioner);

      // Repartitioner
      RCP<RepartitionFactory> repartFactory = rcp(new RepartitionFactory());
      ParameterList repartParams = *(repartFactory->GetValidParameterList());
      MUELU_TEST_AND_SET_PARAM(repartParams, "startLevel",          paramList, defaultList, "repartition: start level",       int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "startLevel",          paramList, defaultList, "repartition: start level",       int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "minRowsPerProcessor", paramList, defaultList, "repartition: min rows per proc", int);
      MUELU_TEST_AND_SET_PARAM(repartParams, "nonzeroImbalance",    paramList, defaultList, "repartition: max imbalance",     double);
      MUELU_TEST_AND_SET_PARAM(repartParams, "remapPartitions",     paramList, defaultList, "repartition: remap parts",       bool);
      MUELU_TEST_AND_SET_PARAM(repartParams, "alwaysKeepProc0",     paramList, defaultList, "repartition: keep proc 0",       bool);
      repartFactory->SetParameterList(repartParams);
      repartFactory->SetFactory("A",         manager.GetFactory("A"));
      repartFactory->SetFactory("Partition", manager.GetFactory("Partition"));
      manager.SetFactory("Importer", repartFactory);

      // Rebalanced A
      RCP<RebalanceAcFactory> newA = rcp(new RebalanceAcFactory());
      newA->  SetFactory("A",         manager.GetFactory("A"));
      newA->  SetFactory("Importer",  manager.GetFactory("Importer"));
      manager.SetFactory("A",         newA);

      // Rebalanced P
      RCP<RebalanceTransferFactory> newP = rcp(new RebalanceTransferFactory());
      ParameterList newPparams;
      newPparams.set("type",     "Interpolation");
      newPparams.set("implicit", this->implicitPRrebalance_);
      newP->  SetParameterList(newPparams);
      newP->  SetFactory("Importer",    manager.GetFactory("Importer"));
      newP->  SetFactory("P",           manager.GetFactory("P"));
      manager.SetFactory("P",           newP);

      // Rebalanced R
      RCP<RebalanceTransferFactory> newR = rcp(new RebalanceTransferFactory());
      ParameterList newRparams;
      newRparams.set("type",     "Restriction");
      newRparams.set("implicit", this->implicitPRrebalance_);
      newR->  SetParameterList(newRparams);
      newR->  SetFactory("Importer",    manager.GetFactory("Importer"));
      newR->  SetFactory("R",           manager.GetFactory("R"));
      newR->  SetFactory("Nullspace",   manager.GetFactory("Ptent"));
      newR->  SetFactory("Coordinates", manager.GetFactory("Coordinates"));
      manager.SetFactory("R",           newR);
      manager.SetFactory("Coordinates", newR);

      // NOTE: the role of NullspaceFactory is to provide nullspace on the finest
      // level if a user does not do that. For all other levels it simply passes
      // nullspace from a real factory to whoever needs it. If we don't use
      // repartitioning, that factory is "TentativePFactory"; if we do, it is
      // "RebalanceTransferFactory". But we still have to have NullspaceFactory as
      // the "Nullspace" of the manager
      nullSpace->SetFactory("Nullspace", newR);
#else
      throw Exceptions::RuntimeError("No repartitioning available for a serial run");
#endif
    }
  }
#undef MUELU_READ_2LIST_PARAM
#undef MUELU_TEST_AND_SET_PARAM

  static bool compare(const ParameterList& list1, const ParameterList& list2) {
    // First loop through and validate the parameters at this level.
    // In addition, we generate a list of sublists that we will search next
    for (ParameterList::ConstIterator it = list1.begin(); it != list1.end(); it++) {
      const std::string&             name   = it->first;
      const Teuchos::ParameterEntry& entry1 = it->second;

      const Teuchos::ParameterEntry *entry2 = list2.getEntryPtr(name);
      if (!entry2)                                           // entry is not present in the second list
        return false;
      if (entry1.isList() && entry2->isList()) {             // sublist check
        compare(Teuchos::getValue<ParameterList>(entry1), Teuchos::getValue<ParameterList>(*entry2));
        continue;
      }
      if (entry1.getAny(false) != entry2->getAny(false))     // entries have different types or different values
        return false;
    }

    return true;
  }
  static bool areSame(const ParameterList& list1, const ParameterList& list2) {
    return compare(list1, list2) && compare(list2, list1);
  }

} // namespace MueLu

#define MUELU_EASYPARAMETERLISTINTERPRETER_SHORT
#endif /* MUELU_EASYPARAMETERLISTINTERPRETER_DEF_HPP */
