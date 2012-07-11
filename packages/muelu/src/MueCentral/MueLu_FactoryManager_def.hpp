#ifndef MUELU_FACTORYMANAGER_DEF_HPP
#define MUELU_FACTORYMANAGER_DEF_HPP

#include "MueLu_FactoryManager_decl.hpp"

#include <Teuchos_ParameterList.hpp>

// Headers for factories used by default:
#include "MueLu_NoFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
//#include "MueLu_GaussSeidelSmoother.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_CoalesceDropFactory2.hpp" // FIXME and default factories
#include "MueLu_RepartitionFactory.hpp"
#include "MueLu_ZoltanInterface.hpp"
#include "MueLu_MultiVectorTransferFactory.hpp"
 
namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::FactoryManager(const RCP<const FactoryBase> PFact, const RCP<const FactoryBase> RFact, const RCP<const FactoryBase> AcFact) 
  { 
    if (PFact  != Teuchos::null) SetFactory("P", PFact);
    if (RFact  != Teuchos::null) SetFactory("R", RFact);
    if (AcFact != Teuchos::null) SetFactory("A", AcFact);

    SetIgnoreUserData(false); // set IgnorUserData flag to false (default behaviour)     
  }
    
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~FactoryManager() { }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetFactory(const std::string & varName, const RCP<const FactoryBase> & factory) {
    if (IsAvailable(varName, factoryTable_)) // TODO: too much warnings (for smoothers)
      GetOStream(Warnings1, 0) << "Warning: FactoryManager::SetFactory(): Changing an already defined factory for " << varName << std::endl;

    factoryTable_[varName] = factory;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const FactoryBase> & FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetFactory(const std::string & varName) const {
    if (FactoryManager::IsAvailable(varName, factoryTable_))
      return factoryTable_.find(varName)->second; // == factoryTable_[varName] (operator std::map[] is not const)
    else 
      return GetDefaultFactory(varName);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const FactoryBase> & FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDefaultFactory(const std::string & varName) const {
    if (IsAvailable(varName, defaultFactoryTable_)) {
        
      return defaultFactoryTable_[varName];
        
    } else {
          
      if (varName == "A")             return SetAndReturnDefaultFactory(varName, rcp(new RAPFactory()));
      if (varName == "P")             return SetAndReturnDefaultFactory(varName, rcp(new SaPFactory(GetFactory("Ptent")))); // GetFactory("Ptent"): Use the same factory instance for both "P" and "Nullspace"
      if (varName == "R")             return SetAndReturnDefaultFactory(varName, rcp(new TransPFactory()));
#if defined(HAVE_MUELU_ZOLTAN) && defined(HAVE_MPI)
//       if (varName == "Partition")     {
//         return SetAndReturnDefaultFactory(varName, rcp(new ZoltanInterface()));
//       }
#endif //ifdef HAVE_MPI
      if (varName == "Importer") {
        return SetAndReturnDefaultFactory(varName, NoFactory::getRCP()); // allows calls to IsAvailable() from RAPFActory
      }
      //JJH FIXME is this going to bite me in the backside?
//       if (varName == "Coordinates") {
//         return SetAndReturnDefaultFactory(varName, rcp(new MueLu::MultiVectorTransferFactory<SC,LO,GO,NO,LMO>(varName,"R")));
//       }

      if (varName == "Nullspace") {
        return SetAndReturnDefaultFactory(varName, rcp(new NullspaceFactory(Teuchos::null, GetFactory("Ptent")))); // GetFactory("Ptent"): Use the same factory instance for both "P" and "Nullspace"
      }

      if (varName == "Graph")               return SetAndReturnDefaultFactory(varName, rcp(new CoalesceDropFactory()));
      if (varName == "UnAmalgamationInfo")  return SetAndReturnDefaultFactory(varName, rcp(new CoalesceDropFactory2()));
      if (varName == "Aggregates")          return SetAndReturnDefaultFactory(varName, rcp(new UCAggregationFactory()));

      // Same factory for both Pre and Post Smoother. Factory for key "Smoother" can be set by users.
      if (varName == "PreSmoother")         return GetFactory("Smoother");
      if (varName == "PostSmoother")        return GetFactory("Smoother");
      
      //if (varName == "Smoother")    return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(rcp(new GaussSeidelSmoother()))));
      if (varName == "Smoother") {
        Teuchos::ParameterList smootherParamList;
        smootherParamList.set("relaxation: type", "Symmetric Gauss-Seidel");
        smootherParamList.set("relaxation: sweeps", (LO) 1);
        smootherParamList.set("relaxation: damping factor", (double) 1.0); //FIXME once Ifpack2's parameter list validator is fixed, change this back to Scalar
        return SetAndReturnDefaultFactory(varName, rcp( new SmootherFactory(rcp(new TrilinosSmoother("RELAXATION", smootherParamList)))));
      }

      if (varName == "CoarseSolver")  return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(rcp(new DirectSolver()),Teuchos::null)));

      if (varName == "Ptent") return SetAndReturnDefaultFactory(varName, rcp(new TentativePFactory())); // Use the same factory instance for both "P" and "Nullspace"

      TEUCHOS_TEST_FOR_EXCEPTION(true, MueLu::Exceptions::RuntimeError, "MueLu::FactoryManager::GetDefaultFactory(): No default factory available for building '"+varName+"'.");
    }
      
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Clean() const { defaultFactoryTable_.clear(); }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const FactoryBase> & FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetAndReturnDefaultFactory(const std::string & varName, const RCP<const FactoryBase> & factory) const {
    TEUCHOS_TEST_FOR_EXCEPTION(factory == Teuchos::null, Exceptions::RuntimeError, "");

    GetOStream(Warnings0,  0) << "Warning: No factory has been specified for building '" << varName << "'." << std::endl;
    GetOStream(Warnings00, 0) << "         Using default factory ";
    { Teuchos::OSTab tab(getOStream(), 7); factory->describe(GetOStream(Warnings00), GetVerbLevel());}

    defaultFactoryTable_[varName] = factory;

    return defaultFactoryTable_[varName];
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  bool FactoryManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::IsAvailable(const std::string & varName, const std::map<std::string, RCP<const FactoryBase> > & factoryTable) {
    return factoryTable.find(varName) != factoryTable.end();
  }

} // namespace MueLu

//TODO: add operator[]
//TODO: should we use a parameterList instead of a std::map? It might be useful to tag which factory have been used and report unused factory.
//TODO: add an option 'NoDefault' to check if we are using any default factory.
//TODO: use Teuchos::ConstNonConstObjectContainer to allow user to modify factories after a GetFactory()

#endif // MUELU_FACTORYMANAGER_DEF_HPP
