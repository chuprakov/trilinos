#ifndef TSFUNSHAREDPTR_H
#define TSFUNSHAREDPTR_H

#include <iostream>

#include "TSFConfig.h"
#include <string>

namespace TSF
{
  using std::string;


  template <class T> class TSFCopier
    {
    public:
      static T* copy(const T& t) {return new T(t);}
    };

  template <class  T> class TSFCloner
    {
    public:
      static T* copy(const T& t) {return t.clone();}
    };




  /**
   * \ingroup Utilities
   * Templated unshared pointer class. Copying an unshared pointer causes a deep copy
   * of the contents.
   */

  template <class T, class C = TSFCopier<T> >
    class TSFUnsharedPtr
    {
    public:
      /** Construct from a shared ptr */
      TSFUnsharedPtr(T* ptr = 0);

      /** Copy ctor, does a deep copy */
      TSFUnsharedPtr(const TSFUnsharedPtr<T, C>& other);

      /** Destructor */
      ~TSFUnsharedPtr();

      /** Assignment op, does a deep copy */
      const TSFUnsharedPtr<T, C>& operator =(const TSFUnsharedPtr<T, C>& rhs);

      /** */
      bool isNull() const;

      /** Access to the raw pointer */
      T* operator ->() const;

      /** Access to the pointer contents */
      T& operator *() const;

      /** Access to the raw pointer */
      T* get() const;

      /** Access to the raw pointer */
      T* ptr() const { return ptr_; }

      /** Report an error. We bypass the ordinary TSF error handling mechanism since it relies
       * on smart pointers */
      void error(const string& msg) const ;

    protected:

      T *ptr_;
    };

  void unsharedPtrError(const string& msg);

  template <class T, class C>
    TSFUnsharedPtr<T, C>::TSFUnsharedPtr(T* ptr)
    : ptr_(ptr)
    {}

  template <class T, class C>
    TSFUnsharedPtr<T, C>::TSFUnsharedPtr(const TSFUnsharedPtr<T, C>& other)
    : ptr_(0)
    {
      if (other.ptr_!=0)
        {
          ptr_ = C::copy(*(other.ptr_));
          if (ptr_==0) error("TSFUnsharedPtr<T, C>::TSFUnsharedPtr copy ctor failed");
        }
    }

  template <class T, class C>
    TSFUnsharedPtr<T, C>::~TSFUnsharedPtr()
    {
      delete ptr_;
    }

  template <class T, class C>
    const TSFUnsharedPtr<T, C>& TSFUnsharedPtr<T, C>::operator =(const TSFUnsharedPtr<T, C>& other)
  {
    if (ptr_ != other.ptr_)
      {
        if (ptr_ != 0) delete ptr_;
        if (other.ptr_!=0)
        {
          ptr_ = C::copy(*(other.ptr_));
          if (ptr_==0) error("TSFUnsharedPtr<T, C>::TSFUnsharedPtr copy ctor failed");
        }
      }
    return *this;
  }

  template <class T, class C>
    bool TSFUnsharedPtr<T, C>::isNull() const {
    return (ptr_ == 0);
  }

  template <class T, class C>
    T* TSFUnsharedPtr<T, C>::operator ->() const {

    if (ptr_ == 0)
      error("TSFUnsharedPtr<T, C>::operator ->() on null pointer");
    return ptr_;
  }

  template <class T, class C>
    T& TSFUnsharedPtr<T, C>::operator *() const {
    if (ptr_ == 0)
      error("TSFUnsharedPtr<T, C>::operator *() on null pointer");
    return *ptr_;
  }

  template <class T, class C>
    T* TSFUnsharedPtr<T, C>::get() const
    {
      return ptr_;
    }

  template <class T, class C>
    void TSFUnsharedPtr<T, C>::error(const string& msg) const
    {
      unsharedPtrError(msg);
    }

  template <class T, class C>
    string toString(const TSFUnsharedPtr<T, C>& sp)
    {
      return "UnsharedPtr<" + TSF::toString(*sp) + ">";
    }
}

#endif
