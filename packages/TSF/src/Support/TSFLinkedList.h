#ifndef LINKLIST_H
#define LINKLIST_H

#include "TSFDefs.h"
#include <iostream>
#include "TSFArray.h"


namespace TSF
{


  /**
   * \ingroup Containers
   * Node in a TSFLinkedList object; for low-level use only.
   */

  template<class T>
    class ListNode {
    public:
    ListNode();
    inline ListNode(const T& data, ListNode<T>* prev);
    inline ListNode(const T& data, ListNode<T>* prev, ListNode<T>* next);
    inline ~ListNode();

    inline ListNode<T>* next() const {return next_;}
    inline ListNode<T>* prev() const {return prev_;}

    inline bool isMaster() const {return isMaster_;}

    inline void setPrev(ListNode<T>* prev) {prev_ = prev;}
    inline void setNext(ListNode<T>* next) {next_ = next;}
    inline void makeMaster() {isMaster_ = true;}

    const T& data() const {return data_;}

    void debug(ostream& os) ;
    private:
    // declare copy ctor and assignment private so they won't compile.
    // We have support for deep copy, but it's handled at the TSFLinkedList level not here.
    ListNode(const ListNode<T>& other);
    ListNode<T>& operator=(const ListNode<T>& other);

    // data members.
    ListNode<T>* next_;
    ListNode<T>* prev_;
    T data_;
    bool isMaster_;
  };

  template<class T> ListNode<T>::ListNode()
    : next_(0), prev_(0), isMaster_(false)
    {;}

  template<class T> inline ListNode<T>::ListNode(const T& data,
                                                 ListNode<T>* prev, ListNode<T>* next)
    : next_(next), prev_(prev), data_(data), isMaster_(false)
    {;}

  template<class T> inline ListNode<T>::ListNode(const T& data, ListNode<T>* prev)
    : next_(0), prev_(prev), data_(data), isMaster_(false)
    {;}


  // recursively delete until we circle back to a master link.

  template<class T> inline ListNode<T>::~ListNode()
    {
      if (next_ != 0 && !next_->isMaster()) delete next_;
    }

  template<class T> inline void ListNode<T>::debug(ostream& os)
    {
      if (isMaster_) cerr << "Master";
      os << "ListNode[ addr = " << this << " p = "
         << prev_ << " n = " << next_ << " data = "
         << data_ << " ]";
    }

  template<class T> ostream& operator<<(ostream& os, ListNode<T>& n)
    {
      n.debug(os);
      return os;
    }

  /**
   * \ingroup Containers
   * Templated linked-list class.
   * Supports prepend and append. Insertion is not needed for now.
   * The data type must have a good empty ctor, copy ctor, and assignment op.
   */

  template<class T>
    class TSFLinkedList {
    public:
    TSFLinkedList();
    inline ~TSFLinkedList();
    TSFLinkedList(const TSFLinkedList<T>& other);
    TSFLinkedList<T>& operator=(const TSFLinkedList<T>& other);

    void prepend(const T& data);
    void append(const T& data);

    bool getFirstNode(T& data)  ;
    bool getLastNode(T& data) ;
    bool getNextNode(T& data)  ;
    bool getPrevNode(T& data) ;

    void removeFirst();
    void removeLast();

    int length() const ;
    private:
    ListNode<T>* master_;
    ListNode<T>* current_;
    int count_;
  };



  template<class T> TSFLinkedList<T>::TSFLinkedList()
    : master_(0), current_(0), count_(0)
    {
      master_ = new ListNode<T>();
      if (master_==0) TSFError::raise("TSFLinkedList empty ctor: can't allocate master");
      master_->makeMaster();
      master_->setPrev(master_);
      master_->setNext(master_);
    }


  template<class T> inline TSFLinkedList<T>::~TSFLinkedList()
    {
      if (master_ != 0) delete master_;
    }


  template<class T> TSFLinkedList<T>::TSFLinkedList(const TSFLinkedList& other)
    : master_(0), current_(0), count_(other.count_)
    {
      master_ = new ListNode<T>();
      if (master_==0) TSFError::raise("TSFLinkedList copy ctor: can't allocate master");
      master_->makeMaster();

      if (other.master_ != 0)
        {
          ListNode<T>* newCurrent = master_;
          ListNode<T>* oldCurrent = other.master_;
          do
            {
              if (oldCurrent->next()->isMaster()) break;
              ListNode<T>* tmp = new ListNode<T>(oldCurrent->next()->data(), newCurrent);
              if (tmp==0) TSFError::raise("TSFLinkedList copy ctor: can't allocate node");
              newCurrent->setNext(tmp);
              newCurrent = newCurrent->next();
              oldCurrent = oldCurrent->next();
            }
          while (!oldCurrent->isMaster());
          // we've looped over all elements. Now close the circle by setting the master's
          // prev pointer to the last element created.
          newCurrent->setNext(master_);
          master_->setPrev(newCurrent);
        }
    }

  template<class T> TSFLinkedList<T>& TSFLinkedList<T>::operator=(const TSFLinkedList<T>& other)
    {
      if (other.master_ == master_) return *this;

      current_ = 0;
      if (master_==0) delete master_;
      master_ = new ListNode<T>();
      if (master_==0) TSFError::raise("TSFLinkedList operator=: can't allocate master");
      master_->makeMaster();

      if (other.master_ != 0)
        {
          ListNode<T>* newCurrent = master_;
          ListNode<T>* oldCurrent = other.master_;
          do
            {
              if (oldCurrent->next()->isMaster()) break;
              ListNode<T>* tmp = new ListNode<T>(oldCurrent->next()->data(), newCurrent);
              if (tmp==0) TSFError::raise("TSFLinkedList operator=: can't allocate node");
              newCurrent->setNext(tmp);
              newCurrent = newCurrent->next();
              oldCurrent = oldCurrent->next();
            }
          while (!oldCurrent->isMaster());
          // we've looped over all elements. Now close the circle by setting the master's
          // prev pointer to the last element created.
          master_->setPrev(newCurrent);
        }

      count_ = other.count_;
      return *this;
    }


  template<class T> void TSFLinkedList<T>::prepend(const T& data)
    {
      if (master_==0) TSFError::raise("bad list master detected");
      ListNode<T>* oldBase = master_->next();
      ListNode<T>* newBase = new ListNode<T>(data, master_, oldBase);
      if (newBase==0) TSFError::raise("TSFLinkedList::prepend()");
      master_->setNext(newBase);
      oldBase->setPrev(newBase);
      count_++;
    }

  template<class T> inline void TSFLinkedList<T>::append(const T& data)
    {
      if (master_==0) TSFError::raise("bad list master detected");
      ListNode<T>* oldEnd = master_->prev();
      ListNode<T>* newEnd = new ListNode<T>(data, oldEnd, master_);
      if (newEnd==0) TSFError::raise("TSFLinkedList::append()");
      master_->setPrev(newEnd);
      oldEnd->setNext(newEnd);
      count_++;
    }

  template<class T> void TSFLinkedList<T>::removeFirst()
    {
      if (master_==0) TSFError::raise("TSFLinkedList<T>::removeFirst() bad master detected");
      if (count_==0) TSFError::raise("TSFLinkedList::removeFirst() called on empty list");
      ListNode<T>* oldBase = master_->next();
      if (current_ == oldBase) current_ = oldBase->next();
      // link up master and the second entry
      master_->setNext(oldBase->next());
      (oldBase->next())->setPrev(master_);
      // set next to zero in removed node, to prevent deletion of whole list.
      oldBase->setNext(0);
      count_--;
      delete oldBase;
    }

  template<class T> void TSFLinkedList<T>::removeLast()
    {
      if (master_==0) TSFError::raise("TSFLinkedList<T>::removeLast() bad master detected");
      if (count_==0) TSFError::raise("TSFLinkedList::removeLast() called on empty list");
      ListNode<T>* oldEnd = master_->prev();
      if (current_ == oldEnd) current_ = oldEnd->prev();
      // link up master and the next-to-last entry
      master_->setPrev(oldEnd->prev());
      (oldEnd->prev())->setNext(master_);
      // set next to zero in removed node, to prevent deletion of whole list.
      oldEnd->setNext(0);
      count_--;
      delete oldEnd;
    }


  template<class T> bool TSFLinkedList<T>::getFirstNode(T& data)
    {
      if (master_==0) TSFError::raise("TSFLinkedList<T>::getFirstNode bad list master detected");
      current_ = master_->next();
      if (current_==0) TSFError::raise("TSFLinkedList<T>::getFirstNode bad list pointer detected");
      if (current_->isMaster()) return false;
      data = current_->data();
      return true;
    }

  template<class T> bool TSFLinkedList<T>::getLastNode(T& data)
    {
      if (master_==0) TSFError::raise("TSFLinkedList<T>::getLastNode bad list master detected");
      current_ = master_->prev();
      if (current_==0) TSFError::raise("TSFLinkedList<T>::getLastNode bad list pointer detected");
      if (current_->isMaster()) return false;
      data = current_->data();
      return true;
    }

  template<class T> bool TSFLinkedList<T>::getNextNode(T& data)
    {
      if (current_==0) TSFError::raise(" TSFLinkedList<T>::getNextNode bad list pointer detected");
      if (current_->next()==0) TSFError::raise(" TSFLinkedList<T>::getNextNode bad list pointer detected");
      if (current_->next()->isMaster()) return false;

      current_ = current_->next();
      data = current_->data();
      return true;
    }

  template<class T> bool TSFLinkedList<T>::getPrevNode(T& data)
    {
      if (current_==0) TSFError::raise(" TSFLinkedList<T>::getNextNode bad list pointer detected");
      if (current_->prev()==0) TSFError::raise(" TSFLinkedList<T>::getNextNode bad list pointer detected");
      if (current_->prev()->isMaster()) return false;

      current_ = current_->prev();
      data = current_->data();
      return true;
    }

  template<class T> inline int TSFLinkedList<T>::length() const
    {
      return count_;
    }


  template<class T> TSFArray<T> arrayify(TSFLinkedList<T>& list)
    {
      TSFArray<T> a(list.length());
      T dummy;

      if (list.length()==0) return a;
      list.getFirstNode(dummy);
      a[0] = dummy;
      int i = 1;
      while (list.getNextNode(dummy))
        {
          a[i] = dummy;
          i++;
        }
      return a;
    }



}
#endif
