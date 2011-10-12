//SKIP
#define MUELU_DEFAULTFACTORYHANDLER_HPP
//SKIP
//SKIP
//SKIP

#ifndef MUELU_DEFAULTFACTORYHANDLER_HPP
#define MUELU_DEFAULTFACTORYHANDLER_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"

#include "MueLu_DefaultFactoryHandlerBase.hpp"

// Headers for factories used by default:
#include "MueLu_SaPFactory.hpp"
//#include "MueLu_TentativePFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_ReUseFactory.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_CoalesceDropFactory.hpp"

namespace MueLu {

  //! Class that provides default factories within Needs class.
  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
  class DefaultFactoryHandler : public DefaultFactoryHandlerBase {
#include "MueLu_UseShortNames.hpp"

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    DefaultFactoryHandler() {}

    //! Destructor.
    virtual ~DefaultFactoryHandler() { }

    //@}
    
    //@{ Get/Set functions.

    virtual const FactoryBase & GetDefaultFactory(const std::string & varName) {
      if (! DefaultFactoryHandlerBase::IsAvailable(varName)) {

        if (varName == "A")            return *NoFactory::get();
        if (varName == "P")            return *NoFactory::get();
        if (varName == "R")            return *NoFactory::get();

        // if (varName == "A")            return SetAndReturnDefaultFactory(varName, rcp(new RAPFactory));
    	if (varName == "Nullspace")    return SetAndReturnDefaultFactory(varName, rcp(new NullspaceFactory()));
        if (varName == "Graph")        return SetAndReturnDefaultFactory(varName, rcp(new CoalesceDropFactory()));
        if (varName == "Aggregates")   return SetAndReturnDefaultFactory(varName, rcp(new UCAggregationFactory()));
//         if (varName == "PreSmoother")  return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(Teuchos::null))); //TODO: change default
//         if (varName == "PostSmoother") return SetAndReturnDefaultFactory(varName, rcp(new SmootherFactory(Teuchos::null)));

        TEST_FOR_EXCEPTION(1, MueLu::Exceptions::RuntimeError, "DefaultFactoryHandler::GetDefaultFactory(): No default factory available for building '"+varName+"'.");
      }

      return DefaultFactoryHandlerBase::GetDefaultFactory(varName);
    }

    //TODO GetDefaultObject, for Smoothers

    //@}
    
  private:

    //! helper
    const FactoryBase & SetAndReturnDefaultFactory(const std::string & varName, const RCP<FactoryBase> factory) {

      GetOStream(Warnings0, 0)  << "Warning: No factory have been specified for building '" << varName << "'." << std::endl;
      GetOStream(Warnings00, 0) << "         using default factory: ";
      { Teuchos::OSTab tab(getOStream(), 8); factory->describe(GetOStream(Warnings00), getVerbLevel()); }

      DefaultFactoryHandlerBase::SetDefaultFactory(varName, factory);
      return DefaultFactoryHandlerBase::GetDefaultFactory(varName); //return factory;
    }

  }; // class DefaultFactoryHandler

} // namespace MueLu

#define MUELU_DEFAULTFACTORYHANDLER_SHORT
#endif //ifndef MUELU_DEFAULTFACTORYHANDLER_HPP

// TODO: print default factory used in debug mode
// TODO: if varName already available in Level but using another factory, print a warning.
// TODO: desallocation of DefaultFactory not set by users
