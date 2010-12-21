#ifndef MUELU_SMOOTHERPROTOTYPE_HPP
#define MUELU_SMOOTHERPROTOTYPE_HPP

#include <iostream>

#include "MueLu_Needs.hpp"
#include "MueLu_Level.hpp"

namespace MueLu {

/*!
  @class SmootherPrototype
  @brief Base class for smoother prototypes

 A smoother prototype is a smoother which can be in two states:
  - ready to be duplicated (parameters defined)
  - ready to be used (setup phase completed)

 'Smoother prototypes' can be fully copied using the Copy() method
 or the copy constructor. They can also copy the parameters of
 another smoother object of the same type (CopyParameters()). Both
 capabilities are used by the SmootherFactory.
*/

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class Level;

template <class ScalarType,class LocalOrdinal,class GlobalOrdinal,class Node, class LocalMatOps>
class SmootherPrototype : public Needs {

#include "MueLu_UseShortNames.hpp"

  private:
    bool IsSetup_;

  public:
    //!@nameConstructors/Destructors.
    //@{
    SmootherPrototype() {}

    virtual ~SmootherPrototype() {}
    //@}

    //! @name Build methods.
    //@{

    virtual void CopyParameters(RCP<SmootherPrototype>) = 0;

    virtual void Setup(RCP<Level<ScalarType, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >) = 0;

    virtual RCP<SmootherPrototype> Copy() = 0;

    //@}

    //! @name Get/Set methods.
    //@{

    //! Get the state of a smoother prototype.
    bool IsSetup() {
      return IsSetup_;
    }

    /*!
      @brief Set the state of a smoother prototype.

      This method should be called by Setup().
    */
    void IsSetup(bool ToF) {
      IsSetup_ = ToF;
    }
    //@}

}; //class SmootherPrototype

} //namespace MueLu

#define MUELU_SMOOTHERPROTOTYPE_SHORT

#endif //ifndef MUELU_SMOOTHERPROTOTYPE_HPP
