#ifndef TSFNONDUPARRAY_H
#define TSFNONDUPARRAY_H

#include "TSFDefs.h"
#include "TSFArray.h"

namespace TSF
{
  using std::string;

  /**
   * \ingroup Containers
   * TSFArray class in which elements cannot be repeated.
   */


  template<class T>
    class TSFNonDupArray
    {
    public:
      /** Create an empty array */
      TSFNonDupArray() : a_(0) {;}
      /** Create from an existing ordinary array, removing any
          duplicates found  */
      TSFNonDupArray(const TSFArray<T>& array);

      /** Preallocate space for n entries */
      void reserve(int n) {a_.reserve(n);}
      /** Get number of elements in the array */
      int length() const {return a_.length();}

      /** check for existence of an entry */
      bool contains(const T& entry, int& index) const ;

      /** Add a new entry, ignoring duplicates */
      TSFNonDupArray<T>& append(const T& entry);
      /** Merge in the elements of an TSFArray, ignoring duplicates */
      TSFNonDupArray<T>& append(const TSFArray<T>& entry);

      // simple accessors.

      /** Read/write access */
      T& operator[](int i) {return a_[i];}
      /** Read-only access */
      const T& operator[](int i) const {return a_[i];}

      /** Write as a String */
      string toString() const {return a_.toString();}

      /** access data as an ordinary array */
      const TSFArray<T>& array() const {return a_;}
    private:
      TSFArray<T> a_;
    };

  /** \relates TSFNonDupArray write to stream */
  template<class T> ostream& operator<<(ostream& os,
                                        const TSFNonDupArray<T>& array);


  template<class T>
    TSFNonDupArray<T>::TSFNonDupArray(const TSFArray<T>& array)
    : a_(0)
    {
      append(array);
    }


  template<class T>
    TSFNonDupArray<T>& TSFNonDupArray<T>::append(const T& entry)
    {
      int index;
      if (!contains(entry, index)) a_.append(entry);
      return *this;
    }

  template<class T>
    TSFNonDupArray<T>& TSFNonDupArray<T>::append(const TSFArray<T>& array)
    {
      reserve(length() + array.length());
      for (int i=0; i<array.length(); i++) append(array[i]);
      return *this;
    }



  template<class T>
    bool TSFNonDupArray<T>::contains(const T& entry, int& index) const
    {
      index = -1;
      for (int i=0; i<length(); i++)
        {
          if (entry==a_[i]) {index = i; return true;}
        }
      return false;
    }

  template<class T> inline string toString(const TSFNonDupArray<T>& a)
    {
      return a.toString();
    }

  template<class T> inline ostream& operator<<(ostream& os,
                                               const TSFNonDupArray<T>& array)
    {
      return os << toString(array);
    }
}
#endif

