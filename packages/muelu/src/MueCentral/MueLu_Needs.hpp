#ifndef MUELU_NEEDS_HPP
#define MUELU_NEEDS_HPP

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_ParameterList.hpp"

#define MueLu_cout(minimumVerbLevel) \
    if (this->getVerbLevel() >= minimumVerbLevel) *(this->getOStream())

namespace MueLu {
  
/*!
 @class Needs
 @brief Class that allows cross-factory communication of data needs.

 Maintains a list of 'Needs' for a given Level. For example, a restriction factory that
 transposes the tentative prolongator 'Needs' the prolongator factory to save this.

 Derives from Teuchos::VerboseObject.

*/
class Needs : public Teuchos::VerboseObject<Needs> {

/*
  //! Friendly print.  FIXME pretty print doesn't work ;(
  friend std::ostream& operator<<(std::ostream &os, Needs const &foo);
*/

  private:

    //! Prior output level
    Teuchos::EVerbosityLevel priorOutputLevel_;
    //! Stores number of outstanding requests for a need.
    Teuchos::ParameterList countTable_;
    //! Stores data associated with a need.
    Teuchos::ParameterList dataTable_;

  protected:
    Teuchos::RCP<Teuchos::FancyOStream> out_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Default constructor.
    Needs() : out_(this->getOStream()) {
      countTable_.setName("countTable");
      dataTable_.setName("dataTable");
    }

    virtual ~Needs() {}
    //@}

    //! @name Access methods
    //@{
    //! Store need label and its associated data. This does not increment the storage counter.
    template <typename T>
    void Save(const std::string ename, const T &entry) {
      if ( !countTable_.isParameter(ename) ) {
        countTable_.set(ename,0);
      }
      dataTable_.set(ename,entry);
    } //Save

    //! Indicate that an object is needed. This increments the storage counter.
    void Request(const std::string ename) {
      int currentCount = countTable_.get(ename,0);
      countTable_.set(ename,++currentCount);
    } //Request


    //! Get data and decrement the storage counter associated with it.
    template <typename T>
    void CheckOut(const std::string ename, T &value) {
      if (!countTable_.isParameter(ename)) {
        std::string msg =  "Checkout: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        value = dataTable_.get<T>(ename);
        int currentCount = countTable_.get(ename,0);
        if (currentCount == 1) {
          countTable_.remove(ename);
          dataTable_.remove(ename);
        } else {
          countTable_.set(ename,--currentCount);
        }
      }
    } //CheckOut

    /*! @brief Get data without decrementing associated storage counter (i.e., read-only access).
    
        This should be used only for local needs or for debugging purposes.
    */
    template <typename T>
    void Examine(const std::string ename, T &value) {
      if (!countTable_.isParameter(ename)) {
        Teuchos::OSTab tab(out_);
        std::string msg =  "Examine: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        value = dataTable_.get<T>(ename);
      }
    } //Examine
    //@}

    //! @name Utilities.
    //@{
    //! NeedsTable raw print
    //FIXME not covered by unit test right now
    void RawPrint(std::ostream &os, Needs const &foo) {
      std::cout << "name | #requests" << std::endl;
      std::cout << "=================" << std::endl;
      std::cout << foo.countTable_ << std::endl << std::endl;
      std::cout << "name | value" << std::endl;
      std::cout << "============" << std::endl;
      std::cout << foo.dataTable_ << std::endl;
    } //RawPrint

    //! Test whether a need's value has been saved.
    bool IsSaved(const std::string ename) {
      if (dataTable_.isParameter(ename)) return true;
      else                               return false;
    }

    //! Test whether a need has been requested.  Note: this tells nothing about whether the need's value exists.
    bool IsRequested(const std::string ename) {
      if (countTable_.isParameter(ename)) return true;
      else                                return false;
    }

    /*! @brief Return the number of outstanding requests for a need.

        Throws a <tt>std::logic_error</tt> exception if the need either hasn't been requested or
        hasn't been saved.
    */
    int NumRequests(const std::string ename) {
      //FIXME should we return 0 instead of throwing an exception?
      if (!countTable_.isParameter(ename)) {
        std::string msg =  "NumRequests: " + ename + " not found in countTable_";
        throw(std::logic_error(msg));
      } else {
        return countTable_.get(ename,0);
      }
    } //NumRequests

    //@}

}; //class Needs

/*
//! NeedsTable pretty print
std::ostream& operator<<(std::ostream &os, Needs const &foo) {
  std::cout << "name  |  #requests  |  value" << std::endl;
  std::cout << "============================" << std::endl;
  for (Teuchos::ParameterList::ConstIterator param=foo.countTable_.begin(); param!=foo.countTable_.end() ; param++) {
    const std::string pname = foo.countTable_.name(param);
    const int& numRequests = foo.countTable_.get<int>(pname);
    const Teuchos::ParameterEntry &someEntry = foo.dataTable_.getEntry(pname);
    std::cout << pname << " | " << numRequests << " | " << someEntry << std::endl;
  }
  return os;
}
*/

} //namespace MueLu

#define MUELU_NEEDS_SHORT

#endif //ifndef MUELU_NEEDS_HPP
