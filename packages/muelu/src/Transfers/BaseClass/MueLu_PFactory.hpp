#ifndef MUELU_PFACTORY_HPP
#define MUELU_PFACTORY_HPP

#include "MueLu_Exceptions.hpp"
#include "MueLu_Level.hpp"

#include "Teuchos_VerboseObject.hpp"
#define MueLu_cout(minimumVerbLevel) \
    if (this->getVerbLevel() >= minimumVerbLevel) *(this->getOStream())

#include <iostream>

namespace MueLu {

/*!
  @class PFactory
  @brief Factory that provides an interface for a concrete implementation of a prolongation operator.

  For a concrete implementation the user has to overwrite the virtual Build method.
*/

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
//class PFactory public OperatorFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
class PFactory : public Teuchos::VerboseObject<PFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> > {
#include "MueLu_UseShortNames.hpp"

  protected:

     bool reUseGraph_;
     bool reUseAggregates_;
     Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    PFactory() : reUseGraph_(false), reUseAggregates_(false), out_(this->getOStream())
    {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "PFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~PFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Abstract Build method.
    */
    virtual bool BuildP(Level &fineLevel, Level &coarseLevel) const = 0;
    //@}

    //! @name Set methods.
    //@{

    void ReUseAggregates(bool const &value) {
      reUseAggregates_ = value;
    }

    void ReUseGraph(bool const &value) {
      reUseGraph_ = value;
    }
    //@}

    //! @name Get methods.
    //@{

    bool ReUseAggregates() const {
      return reUseAggregates_;
    }

    bool ReUseGraph() const {
      return reUseGraph_;
    }

    //@}

}; //class PFactory

} //namespace MueLu

#define MUELU_PFACTORY_SHORT

#endif //ifndef MUELU_PFACTORY_HPP
