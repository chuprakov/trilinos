#ifndef MUELU_FACTORYFACTORY_DECL_HPP
#define MUELU_FACTORYFACTORY_DECL_HPP

#include <string>

#include <Teuchos_ParameterEntry.hpp>
#include <Teuchos_Array.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_FactoryFactory_fwd.hpp"

#include "MueLu_HierarchyFactory.hpp"

#include "MueLu_FactoryBase.hpp"
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Hierarchy_fwd.hpp"
#include "MueLu_FactoryManagerBase_fwd.hpp"
#include "MueLu_FactoryManager_fwd.hpp"
#include "MueLu_FactoryBase_fwd.hpp"

#include "MueLu_CoalesceDropFactory.hpp" //TMP
#include "MueLu_RAPFactory.hpp" //TMP
#include "MueLu_TransPFactory.hpp" //TMP
#include "MueLu_SaPFactory.hpp" //TMP
#include "MueLu_TrilinosSmoother.hpp" //TMP
#include "MueLu_SmootherFactory.hpp" //TMP
#include "MueLu_TentativePFactory.hpp" //TMP
#include "MueLu_UCAggregationFactory.hpp" //TMP
#include "MueLu_DirectSolver.hpp" //TMP
#include "MueLu_Exceptions.hpp" //TMP
#include "MueLu_MultiVectorTransferFactory.hpp"
#include "MueLu_PermutedTransferFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_RepartitionFactory.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void, LocalOrdinal, Node>::SparseOps>
  class FactoryFactory {
#undef MUELU_FACTORYFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; // TODO: remove

  public:

    // Parameter List Parsing:
    // ---------
    //     <Parameter name="smootherFact0" type="string" value="TrilinosSmoother"/>
    //
    // or:
    //
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       ...
    //     </ParameterList>
    //
    RCP<const FactoryBase> BuildFactory(const Teuchos::ParameterEntry & param, const FactoryMap & factoryMapIn) const {

      //TODO: add test restricted keyword

      // Find factory
      std::string factoryName;
      Teuchos::ParameterList paramList;
      if (!param.isList()) {
        factoryName = Teuchos::getValue<std::string>(param);
      } else {
        paramList = Teuchos::getValue<Teuchos::ParameterList>(param);
        factoryName = paramList.get<std::string>("factory");
      } 

      // TODO: see how Teko handles this.
      if (factoryName == "CoalesceDropFactory") {
        return BuildCoalesceDropFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TentativePFactory") {
        return BuildTentativePFactory(paramList, factoryMapIn);
      }
      if (factoryName == "SaPFactory") {
        return BuildSaPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TransPFactory") {
        return BuildTransPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "RAPFactory") {
        return BuildRAPFactory(paramList, factoryMapIn);
      }
      if (factoryName == "UCAggregationFactory") {
        return BuildUCAggregationFactory(paramList, factoryMapIn);
      }
      if (factoryName == "TrilinosSmoother") {
        return BuildTrilinosSmoother(paramList, factoryMapIn);
      }
      if (factoryName == "DirectSolver") {
        return BuildDirectSolver(paramList, factoryMapIn);
      }
      if (factoryName == "MultiVectorTransferFactory") {
        return BuildMultiVectorTransferFactory(paramList, factoryMapIn);
      }
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
      if (factoryName == "ZoltanInterface") {
        return BuildZoltanInterface(paramList, factoryMapIn);
      }
#endif
      if (factoryName == "RepartitionFactory") {
        return BuildRepartitionFactory(paramList, factoryMapIn);
      }
      if (factoryName == "PermutedTransferFactory") {
        return BuildPermutedTransferFactory(paramList, factoryMapIn);
      }

      // Use a user defined factories (in <Factories> node)
      if (factoryMapIn.find(factoryName) != factoryMapIn.end()) {
        TEUCHOS_TEST_FOR_EXCEPTION((param.isList() && (++paramList.begin() != paramList.end())), Exceptions::RuntimeError, 
                                   "MueLu::FactoryFactory: Error during the parsing of: " << std::endl << paramList << std::endl
                                   << "'" << factoryName << "' is not a factory name but an existing instance of a factory." << std::endl
                                   << "Extra parameters cannot be specified after the creation of the object." << std::endl << std::endl
                                   << "Correct syntaxes includes:" << std::endl
                                   << " <Parameter name=\"...\" type=\"string\" value=\"" << factoryName << "\"/>" << std::endl
                                   << "or" << std::endl
                                   << " <ParameterList name=\"...\"><Parameter name=\"factory\" type=\"string\" value=\"" << factoryName << "\"/></ParameterList>" << std::endl
                                   );

        return factoryMapIn.find(factoryName)->second;
      }

      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory: unknown factory name : " << factoryName);

      return Teuchos::null;
    }
  
    //
    //
    //

    // FOLLOWING FUNCTIONS SHOULD LIVE WITH THE CORRESPONDING CLASS

    //
    //
    //

#define MUELU_FACTORY_PARAM(name, var)                                  \
    RCP<const FactoryBase> var; if (paramList.isParameter(name)) { var = BuildFactory(paramList.getEntry(name), factoryMapIn); }
    
    //! CoalesceDropFactory
    RCP<FactoryBase> BuildCoalesceDropFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      MUELU_FACTORY_PARAM("A", AFact);
      MUELU_FACTORY_PARAM("Nullspace", NullspaceFact);

      return rcp(new CoalesceDropFactory(AFact, NullspaceFact));
    }

    //! TentativePFactory
    RCP<FactoryBase> BuildTentativePFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      MUELU_FACTORY_PARAM("Aggregates", AggFact);
      MUELU_FACTORY_PARAM("Nullspace",  NullspaceFact);
      MUELU_FACTORY_PARAM("A", AFact);

      return rcp(new TentativePFactory(AggFact, NullspaceFact, AFact));
    }
    
    //! SaPFactory
    RCP<FactoryBase> BuildSaPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new SaPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "SaPFactory", Exceptions::RuntimeError, "");      
      MUELU_FACTORY_PARAM("InitialP", InitialPFact);
      MUELU_FACTORY_PARAM("A", AFact);

      RCP<SaPFactory> f = rcp(new SaPFactory(InitialPFact, AFact));
      if (paramList.isParameter("DampingFactor")) f->SetDampingFactor(paramList.get<Scalar>("DampingFactor"));

      return f;
    }

    //! TransPFactory
    RCP<FactoryBase> BuildTransPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
//       if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
//         return rcp(new TransPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TransPFactory", Exceptions::RuntimeError, "");      
      MUELU_FACTORY_PARAM("P", PFact);

      return rcp(new TransPFactory(PFact));
    }

    //! RaPFactory
    RCP<FactoryBase> BuildRAPFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new RAPFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RAPFactory", Exceptions::RuntimeError, "");      
      MUELU_FACTORY_PARAM("P", PFact);
      MUELU_FACTORY_PARAM("R", RFact);
      MUELU_FACTORY_PARAM("A", AFact);

      RCP<RAPFactory> r = rcp(new RAPFactory(PFact, RFact, AFact));

      if (paramList.isSublist("TransferFactories")) {
        Teuchos::ParameterList transferList = paramList.sublist("TransferFactories");
        for (Teuchos::ParameterList::ConstIterator param = transferList.begin(); param != transferList.end(); ++param) {
          RCP<const FactoryBase> p = BuildFactory(transferList.entry(param), factoryMapIn);
          r->AddTransferFactory(p);
        }
      }

      return r;
    }

    //! UCAggregationFactory
    RCP<FactoryBase> BuildUCAggregationFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end()) // short-circuit. Use default parameters of constructor
        return rcp(new UCAggregationFactory());

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "UCAggregationFactory", Exceptions::RuntimeError, "");      
      MUELU_FACTORY_PARAM("Graph", GraphFact);

      RCP<UCAggregationFactory> f = rcp(new UCAggregationFactory(GraphFact));

      if(paramList.isParameter("Ordering")) {
        std::string orderingStr = paramList.get<std::string>("Ordering");
        Ordering ordering;
        if (orderingStr == "Natural")
          ordering = NATURAL;
        else if (orderingStr == "Random")
          ordering = RANDOM;
        else if (orderingStr == "Graph")
          ordering = GRAPH;
        else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::FactoryFactory::BuildUCAggregationFactory()::Unknown Ordering type");

        f->SetOrdering(ordering);
      }

      if(paramList.isParameter("MaxNeighAlreadySelected")) {
        f->SetMaxNeighAlreadySelected(paramList.get<int>("MaxNeighAlreadySelected"));
      }
      
      if(paramList.isParameter("Phase3AggCreation")) {
        f->SetPhase3AggCreation(paramList.get<double>("Phase3AggCreation"));
      }

      if(paramList.isParameter("MinNodesPerAggregate")) {
        f->SetMinNodesPerAggregate(paramList.get<int>("MinNodesPerAggregate"));
      }

      return f;
    }

    //! TrilinosSmoother
    // Parameter List Parsing:
    //     <ParameterList name="smootherFact1">
    //       <Parameter name="factory" type="string" value="TrilinosSmoother"/>
    //       <Parameter name="verbose" type="string" value="Warnings"/>
    //       <Parameter name="type" type="string" value="RELAXATION"/>
    //       <ParameterList name="ParameterList">
    //       ...
    //       </ParameterList>
    //     </ParameterList>
    RCP<FactoryBase> BuildTrilinosSmoother(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new TrilinosSmoother())));

      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "TrilinosSmoother", Exceptions::RuntimeError, "");

      // Is it true? TEUCHOS_TEST_FOR_EXCEPTION(!paramList.isParameter("type"), Exceptions::RuntimeError, "TrilinosSmoother: parameter 'type' is mandatory");
      // type="" is default in TrilinosSmoother, but what happen then?

      std::string type="";            if(paramList.isParameter("type"))          type    = paramList.get<std::string>("type");
      int         overlap=0;          if(paramList.isParameter("overlap"))       overlap = paramList.get<int>        ("overlap");
      // std::string verbose;         if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params;  if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");

      return rcp(new SmootherFactory(rcp(new TrilinosSmoother(type, params, overlap))));
    }
    
    RCP<FactoryBase> BuildDirectSolver(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      if (paramList.begin() == paramList.end())
        return rcp(new SmootherFactory(rcp(new DirectSolver())));
      
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "DirectSolver", Exceptions::RuntimeError, "");

      std::string type;              if(paramList.isParameter("type"))          type = paramList.get<std::string>("type");
      // std::string verbose;        if(paramList.isParameter("verbose"))       verbose = paramList.get<std::string>("verbose");
      Teuchos::ParameterList params; if(paramList.isParameter("ParameterList")) params  = paramList.get<Teuchos::ParameterList>("ParameterList");
      
      return rcp(new SmootherFactory(rcp(new DirectSolver(type, params))));
    }

    RCP<FactoryBase> BuildMultiVectorTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "MultiVectorTransferFactory", Exceptions::RuntimeError, "");      
      
      std::string vectorName="";      vectorName = paramList.get<std::string>("vectorName");
      std::string restrictionName=""; restrictionName = paramList.get<std::string>("restrictionName");
      MUELU_FACTORY_PARAM("R", RFact);
      
      return rcp(new MultiVectorTransferFactory(vectorName, restrictionName, RFact));
    }

#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
    RCP<FactoryBase> BuildZoltanInterface(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "ZoltanInterface", Exceptions::RuntimeError, "");      
      
      MUELU_FACTORY_PARAM("A", AFact);
      MUELU_FACTORY_PARAM("TransferFactory", TransferFactory);
      
      return rcp(new ZoltanInterface(AFact, TransferFactory));
    }
#endif

    RCP<FactoryBase> BuildRepartitionFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "RepartitionFactory", Exceptions::RuntimeError, "");      
      
      MUELU_FACTORY_PARAM("Factory", Factory);
      MUELU_FACTORY_PARAM("A", AFact);
      int minRowsPerProc=2000;     if(paramList.isParameter("minRowsPerProc"))   minRowsPerProc = paramList.get<int>("minRowsPerProc");
      double nonzeroImbalance=1.2; if(paramList.isParameter("nonzeroImbalance")) nonzeroImbalance = paramList.get<double>("nonzeroImbalance");
      int startLevel=1;            if(paramList.isParameter("startLevel"))       startLevel = paramList.get<int>("startLevel");
      int diffusive=0;             if(paramList.isParameter("diffusive"))        diffusive = paramList.get<int>("diffusive");
      
      return rcp(new RepartitionFactory(Factory, AFact, minRowsPerProc, nonzeroImbalance, startLevel, diffusive));
    }

    RCP<FactoryBase> BuildPermutedTransferFactory(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn) const {
      TEUCHOS_TEST_FOR_EXCEPTION(paramList.get<std::string>("factory") != "PermutedTransferFactory", Exceptions::RuntimeError, "");      

      MUELU_FACTORY_PARAM("RepartitionFactory", RepartitionFact);
      MUELU_FACTORY_PARAM("A", AFact);
      MUELU_FACTORY_PARAM("P", PFact);
      
      std::string type; type = paramList.get<std::string>("type");
      if (type == "Interpolation") {
        return rcp(new PermutedTransferFactory(RepartitionFact, AFact, PFact, MueLu::INTERPOLATION));
      } else if (type == "Restriction") {
        MUELU_FACTORY_PARAM("R", RFact);
        MUELU_FACTORY_PARAM("TransferFactory", TransferFactory);
        return rcp(new PermutedTransferFactory(RepartitionFact, AFact, RFact, MueLu::RESTRICTION, PFact, TransferFactory));
      } else {
        TEUCHOS_TEST_FOR_EXCEPT(1);
      }
    }
  }; // class

} // namespace MueLu

#define MUELU_FACTORYFACTORY_SHORT
#endif // MUELU_FACTORYFACTORY_DECL_HPP

// TODO: handle factory parameters
// TODO: parameter validator
// TODO: static
// TODO: default parameters should not be duplicated here and on the Factory (ex: default for overlap (=0) is defined both here and on TrilinosSmoother constructors)
