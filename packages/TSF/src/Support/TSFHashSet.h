#ifndef HASHSET_H
#define HASHSET_H

#include "TSFDefs.h"
#include "TSFArray.h"
#include "TSFHashUtils.h"

namespace TSF
{
  using std::string;


  /** \ingroup Containers
   * TSFHashSet is a hashtable-based set, similar to the STL set class
   * or the Java TSFHashSet class.
   */


  template<class Key> class TSFHashSet
    {
    public:
      /**
       * Create an empty TSFHashSet
       */
      inline TSFHashSet(int capacity=101);

      /**
       * Check for the presence of a key
       */
      inline bool containsKey(const Key& key) const ;

      /**
       * Put a new object into the table.
       */
      inline void put(const Key& key);

      /**
       * Remove from the table the element given by key.
       */
      inline void remove(const Key& key);

      /**
       * Get the number of elements in the table
       */
      inline int size() const {return count_;}

      /**
       * Get list of keys in TSFArray form
       */
      inline TSFArray<Key> arrayify() const ;

      /**
       * render as a string
       */
      string toString() const ;
    private:
      /** rebuild the hashtable when the size has changed */
      inline void rehash();
      /** get the next prime number near a given capacity */
      inline int nextPrime(int newCap) const ;

      TSFArray<TSFArray<Key> > data_;
      int count_;
      int capacity_;
      mutable Key mostRecentKey_;
    };


  /** \relates TSFHashSet write to a stream */
  template<class Key>
    ostream& operator<<(ostream& os, const TSFHashSet<Key>& h);

  template<class Key> inline
    string toString(const TSFHashSet<Key>& h) {return h.toString();}


  template<class Key> inline
    TSFHashSet<Key>::TSFHashSet(int capacity)
    : data_(), count_(0), capacity_(TSFHashUtils::nextPrime(capacity))
    {
      data_.resize(capacity_);
    }

  template<class Key> inline
    bool TSFHashSet<Key>::containsKey(const Key& key) const
    {
      const TSFArray<Key>& candidates
        = data_[hashCode(key) % capacity_];

      for (int i=0; i<candidates.length(); i++)
        {
          const Key& c = candidates[i];
          if (c == key)
            {
              return true;
            }
        }
      return false;
    }

  template<class Key> inline
    void TSFHashSet<Key>::put(const Key& key)
    {
      int index = hashCode(key) % capacity_;

      TSFArray<Key>& local = data_[index];

      // check for duplicate key
      for (int i=0; i<local.length(); i++)
        {
          if (local[i] == key)
            {
              return;
            }
        }

      // no duplicate key, so increment element count by one.
      count_++;

      // check for need to resize.
      if (count_ > capacity_)
        {
          capacity_ = TSFHashUtils::nextPrime(capacity_+1);
          rehash();
          // recaluate index
          index = hashCode(key) % capacity_;
        }

      data_[index].append(key);
    }



  template<class Key> inline
    void TSFHashSet<Key>::rehash()
    {
      TSFArray<TSFArray<Key> > tmp(capacity_);

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              int newIndex = hashCode(data_[i][j]) % capacity_;
              tmp[newIndex].append(data_[i][j]);
            }
        }

      data_ = tmp;
    }

  template<class Key> inline
    TSFArray<Key> TSFHashSet<Key>::arrayify() const
    {
      TSFArray<Key> rtn;
      rtn.reserve(size());

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              rtn.append(data_[i][j]);
            }
        }

      return rtn;
    }

  template<class Key>  inline
    string TSFHashSet<Key>::toString() const
    {
      string rtn = "TSFHashSet[";

      bool first = true;

      for (int i=0; i<data_.length(); i++)
        {
          for (int j=0; j<data_[i].length(); j++)
            {
              if (!first) rtn += ", ";
              first = false;
              rtn += toString(data_[i][j]);
            }
        }
      rtn += "]";
      return rtn;
    }


  template<class Key> inline
    void TSFHashSet<Key>::remove(const Key& key)
    {
      if (!containsKey(key)) TSFError::raise("key not found in hash set");

      count_--;
      int h = hashCode(key) % capacity_;
      TSFArray<Key>& candidates = data_[h];

      for (int i=0; i<candidates.length(); i++)
        {
          if (candidates[i] == key)
            {
              candidates.remove(i);
              break;
            }
        }
    }



  template<class Key>  inline
    ostream& operator<<(ostream& os, const TSFHashSet<Key>& h)
    {
      return os << h.toString();
    }


}

#endif
