#ifndef MUELU_SMOOTHERBASE_HPP
#define MUELU_SMOOTHERBASE_HPP

#include <iostream>

#include "MueLu_Needs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

  using Teuchos::RCP;

/*!
  @class SmootherBase
  @brief Base class for smoothers

  This has the signature for the required Apply function and contains data that is generic across all
  smoothers.
*/

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class SmootherBase : public Needs {

#include "MueLu_UseShortNames.hpp"

  private:
  std::string Type_;

  public:
    //@{ Constructors/Destructors.
    SmootherBase() {}

    virtual ~SmootherBase() {}
    //@}

    //! @name Apply methods.
    //@{

    //! Apply smoother.
    virtual void Apply(RCP<MultiVector> x, RCP<MultiVector> const rhs, bool InitialGuessIsZero) = 0;

    //@}

    //! @name Set/Get methods.
    //@{

    //! Prototype of the method to setup the number of iteration of the smoother.
    virtual void SetNIts(LO Nits) = 0;

    //! Get the smoother type.
    std::string GetType() {
      return Type_;
    }

    /*! @brief Set the smoother type.

      This method must be called by constructors of derived classes.
    */
    void SetType(std::string type) {
      Type_ = type;
    }

    //! @name Utilities.
    //@{

    //! @brief Print information about the smoother
    virtual void Print(std::string prefix) = 0;

    //@}

}; //class SmootherBase

} //namespace MueLu

#define MUELU_SMOOTHERBASE_SHORT

#endif //ifndef MUELU_SMOOTHERBASE_HPP
