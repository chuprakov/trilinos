/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFHANDLE_HPP
#define TSFHANDLE_HPP

#include "TSFConfigDefs.hpp"
#include "TSFPrintable.hpp"
#include "TSFDescribable.hpp"
#include "TSFHandleable.hpp"
#include "Teuchos_RefCountPtr.hpp"

namespace TSFExtended
{
  using namespace Teuchos;

  /**
   * Class TSFExtended::Handle provides a general implementation
   * of the common features of reference-counted handles.
   */
  template <class PointerType>
  class Handle
  {
  public:
    /** Empty ctor  */
    Handle() : ptr_() {;}

    /** Construct */
    Handle(const RefCountPtr<PointerType>& ptr) : ptr_(ptr) {;}

    /** */
    Handle(Handleable<PointerType>* rawPtr) : ptr_(rawPtr->getRcp()) {;}

    /** Read-only access to the underlying smart pointer. */
    const RefCountPtr<PointerType>& ptr() const {return ptr_;}

    /** Read-write access to the underlying smart pointer. */
    RefCountPtr<PointerType> ptr() {return ptr_;}

    /** 
     * Print to a stream using the Printable interface. 
     * If the contents of the handle cannot be 
     * downcasted or crosscasted to a Printable*, an exception
     * will be thrown 
     */
    void print(std::ostream& os) const ;

    /** 
     * Return a short descriptive string using the Describable interface.
     * If the contents of the handle cannot be 
     * downcasted or crosscasted to a Describable*, an exception
     * will be thrown. 
     */
    std::string describe() const ;

  private:
    RefCountPtr<PointerType> ptr_;
  };

  /* implementation of print() */
  template <class PointerType> inline 
  void Handle<PointerType>::print(std::ostream& os) const 
  {
    const Printable* p = dynamic_cast<const Printable*>(ptr_.get());
      
    TEST_FOR_EXCEPTION(p==0, std::runtime_error,
                       "Attempted to cast non-printable "
                       "pointer to a Printable");
      p->print(os);
  }

  /* implementation of describe() */
  template <class PointerType> inline
  std::string Handle<PointerType>::describe() const 
  {
    const Describable* p = dynamic_cast<const Describable*>(ptr_.get());
    
    TEST_FOR_EXCEPTION(p==0, std::runtime_error,
                       "Attempted to cast non-describable "
                       "pointer to a Describable");
    return p->describe();
  }
}


template <class PointerType> inline
ostream& operator<<(ostream& os, const TSFExtended::Handle<PointerType>& h)
{
  h.print(os);
  return os;
}

#endif
