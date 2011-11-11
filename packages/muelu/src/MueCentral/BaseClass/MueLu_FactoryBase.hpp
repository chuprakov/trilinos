#ifndef MUELU_FACTORYBASE_HPP
#define MUELU_FACTORYBASE_HPP

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_Level_fwd.hpp"

namespace MueLu {

  static int generateUniqueFactoryId() {
    static int i = 0;
    ++i;
    return i;
  }
  
  //! Base class for factories (e.g., R, P, and A_coarse).
  class FactoryBase : public BaseClass {

  public:
    //@{ Constructors/Destructors.

    //! Constructor.
    FactoryBase()
      : id_(MueLu::generateUniqueFactoryId())
    { }

    //! Destructor.
    virtual ~FactoryBase() { }
    //@}

    //@{
    //! @name Build methods.

    virtual void CallBuild(Level & requestedLevel) const = 0;

    virtual void CallDeclareInput(Level & requestedLevel) const = 0;
    //@}

    //@{
    //! @name Access factory properties

    /// return unique factory id
    int getID() const { return id_; };

    //@}

  private:
    const int id_;

  }; //class FactoryBase

} //namespace MueLu

#define MUELU_FACTORYBASE_SHORT
#endif //ifndef MUELU_FACTORYBASE_HPP

//TODO: use unique ID instead of ptr in Level
