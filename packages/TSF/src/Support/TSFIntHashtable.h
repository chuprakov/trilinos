#ifndef INTHASHTABLE_H
#define INTHASHTABLE_H

#include "TSFDefs.h"

#include "TSFArray.h"
#include "TSFHashUtils.h"


namespace TSF
{

  class IntPair
    {
    public:
      IntPair() : key_(0), value_(0) {;}
      IntPair(int key, int value) : key_(key), value_(value) {;}

      int key_;
      int value_;
    };

  /**
   * \ingroup Containers
   * TSFHashtable hardwired for integers
   */

  class TSFIntHashtable
    {
    public:
      TSFIntHashtable(int capacity=101);

      bool containsKey(int key) const ;
      int get(int key) const ;

      void put(int key, int value) ;

      int size() const {return count_;}


      string toString() const ;
      const TSFArray<TSFArray<IntPair> >& data() const {return data_;}
    private:
      void rehash();
      int nextPrime(int newCap) const ;

      TSFArray<TSFArray<IntPair> > data_;
      int count_;
      int capacity_;
      mutable int mostRecentValue_;
      mutable int mostRecentKey_;
    };

  inline string toString(const TSFIntHashtable& h)
    {
      return h.toString();
    }


  inline ostream& operator<<(ostream& os, const TSFIntHashtable& h)
    {
      return os << h.toString();
    }

}



#endif
