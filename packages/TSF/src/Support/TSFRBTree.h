#ifndef TSFRBTREE_H
#define TSFRBTREE_H

#include "TSFDefs.h"
#include "TSFRBTreeNode.h"
#include "TSFArray.h"

namespace TSF
{


  /**
   * \ingroup Containers
   * TSFRBTree is an associative array implemented as a red-black binary tree.
   *
   * Red-black trees are described in _Algorithms_ by Cormen, Leiserson, and
   * Rivest. This code closely follows their pseudocode.
   *
   * To use this code on classes Key and Value:
   * (1) Key and Value must have empty ctors
   * (2) Key must have a less-than comparison operator<(const Key& other)
   * (3) To print an TSFRBTree, Key and value should
   * have stream insertion operators.
   *
   */
  template<class Key, class Value> class TSFRBTree
    {
    public:
      /** Create an empty TSFRBTree */
      inline TSFRBTree();
      inline TSFRBTree(const TSFRBTree<Key, Value>& other);
      inline ~TSFRBTree();
      inline const TSFRBTree<Key, Value>& operator=
        (const TSFRBTree<Key, Value>& other);

      /** get number of elements in tree */
      int size() const {return count_;}
      /** get value indexed by key */
      inline const Value& get(const Key& key) const ;
      inline bool get(const Key& key, Value& value) const ;
      /** find out if tree contains an element indexed by key */
      inline bool containsKey(const Key& key) const ;
      /** put a new key,value pair in the tree */
      inline bool put(const Key& key, const Value& value);

      inline ostream& print(ostream& os, bool debug=true) const ;

      /** get arrays of keys and values */
      inline void arrayify(TSFArray<Key>& keys, TSFArray<Value>& values) const ;
      inline TSFArray<Key> arrayifyKeys() const ;
      inline TSFArray<Value> arrayifyValues() const ;
    private:
      inline const TSFRBTreeNode<Key, Value>* minimum() const ;
      inline const TSFRBTreeNode<Key, Value>* maximum() const ;
      inline void rotateLeft(TSFRBTreeNode<Key, Value>* pivot);
      inline void rotateRight(TSFRBTreeNode<Key, Value>* pivot);
      inline const TSFRBTreeNode<Key, Value>* pointerTo(const Key& key) const ;
      inline TSFRBTreeNode<Key, Value>* rawInsert(const Key& key, const Value& value);
      TSFRBTreeNode<Key, Value>* root_;
      int count_;
      Key cacheKey_;
      Value cacheValue_;
    };



  //----------------------- implementation -----------------------------------

  template<class Key, class Value> inline TSFRBTree<Key, Value>::TSFRBTree()
    : root_(0),
    count_(0),
    cacheKey_(),
    cacheValue_()
    {
      ;
    }


  template<class Key, class Value> inline
    TSFRBTree<Key, Value>::TSFRBTree(const TSFRBTree<Key, Value>& other)
    : root_(0),
    count_(0),
    cacheKey_(other.cacheKey_),
    cacheValue_(other.cacheValue_)
    {
      if (other.root_ != 0)
        {
          root_ = other.root_->clone();
          if (root_==0) TSFError::raise("TSFRBTree copy ctor");
          count_ = other.count_;
        }
    }

  template<class Key, class Value> inline TSFRBTree<Key, Value>::~TSFRBTree()
    {
      delete root_;
    }

  template<class Key, class Value> inline
    const TSFRBTree<Key, Value>& TSFRBTree<Key, Value>::operator=
    (const TSFRBTree<Key, Value>& other)
    {
      if (this != &other)
        {
          cacheKey_ = other.cacheKey_;
          cacheValue_ = other.cacheValue_;
          if (root_)
            {
              delete root_;
              root_ = 0;
            }
          if (other.root_ != 0)
            {
              root_ = other.root_->clone();
              if (root_==0) TSFError::raise("TSFRBTree copy ctor");
              count_ = other.count_;
            }
        }
      return *this;
    }

  // -------------------------------------------------------------------------
  //
  //     Insertion: put() and rawInsert()
  //
  //--------------------------------------------------------------------------


  // rawInsert() sticks a node into the tree without regard to red/black
  // coloring. put() does a raw insert, followed by recoloring and rotations
  // for balance.

  template<class Key, class Value> inline
    TSFRBTreeNode<Key, Value>* TSFRBTree<Key, Value>::rawInsert(const Key& key,
                                                                const Value& value)
    {
      count_++;
      // If this is the first entry into the tree, root_ will not exist.
      // In that case, create root with a null parent.
      if (root_==0)
        {
          root_ = new TSFRBTreeNode<Key, Value>(key, value, 0, TSFRBRed);
          return root_;
        }

      // Otherwise, begin at root and do an ordered search for a place to fit
      // the new node.

      TSFRBTreeNode<Key, Value>* newParent = 0 ;
      TSFRBTreeNode<Key, Value>* newPtr = root_;

      // Compare to key at each node along path, and go left/right accordingly.
      // At each turn, keep track of whether the turn was left or right; that
      // will be needed in the final step.
      bool lastWasLeft = false;
      while (newPtr != 0)
        {
          newParent = newPtr;
          if (key < newPtr->key_)
            {
              newPtr = newPtr->left_;
              lastWasLeft = true;
            }
          else if (key == newPtr->key_)
            {
              count_--;
              newPtr->value_ = value;
              return newPtr;
            }
          else
            {
              newPtr = newPtr->right_;
              lastWasLeft = false;
            }
        }

      // OK, we've found a location for the new node.
      // Create a node object, and inform the parent that it has a new
      // left or right child.
      TSFRBTreeNode<Key, Value>* newNode
        = new TSFRBTreeNode<Key, Value>(key, value, newParent, TSFRBRed);
      if (lastWasLeft)
        {
          newParent->left_ = newNode;
        }
      else
        {
          newParent->right_ = newNode;
        }

      // all done; return a ptr to the newly created node.
      return newNode;
    }


  // put(). Do a rawInsert(), then fix up balance.

  template<class Key, class Value> inline
    bool TSFRBTree<Key, Value>::put(const Key& key, const Value& value)
    {
      // stick the new node in, coloring it red.
      TSFRBTreeNode<Key, Value>* x = rawInsert(key, value);
      if (x==0) return false;

      // now fix up the RB state of the tree.

      while (x != root_ && x->parent_ != 0 && x->parent_->color_ == TSFRBRed)
        {
          if (x->grandparent() == 0) TSFError::raise("TSFRBTree::put() null grandparent discovered");
          if (x->parent_->isLeftChild())
            {
              TSFRBTreeNode<Key, Value>* y = x->grandparent()->right_;
              if (y!=0  && y->color_ == TSFRBRed)
                {
                  x->parent_->color_ = TSFRBBlack;
                  y->color_ = TSFRBBlack;
                  x->grandparent()->color_ = TSFRBRed;
                  x = x->grandparent();
                }
              else if (!x->isLeftChild())
                {
                  x = x->parent_;
                  rotateLeft(x);
                }
              else
                {
                  x->parent_->color_ = TSFRBBlack;
                  x->grandparent()->color_ = TSFRBRed;
                  rotateRight(x->grandparent());
                }
            }
          else
            {
              TSFRBTreeNode<Key, Value>* y = x->grandparent()->left_;
              if (y!=0 && y->color_ == TSFRBRed)
                {
                  x->parent_->color_ = TSFRBBlack;
                  y->color_ = TSFRBBlack;
                  x->grandparent()->color_ = TSFRBRed;
                  x = x->grandparent();
                }
              else if (x->isLeftChild())
                {
                  x = x->parent_;
                  rotateRight(x);
                }
              else
                {
                  x->parent_->color_ = TSFRBBlack;
                  x->grandparent()->color_ = TSFRBRed;
                  rotateLeft(x->grandparent());
                }
            }
        }

      root_->color_ = TSFRBBlack;
      return true;
    }






  template<class Key, class Value> inline
    const TSFRBTreeNode<Key, Value>* TSFRBTree<Key, Value>::minimum() const
    {
      if (root_==0) return 0;
      else return root_->minimum();
    }

  template<class Key, class Value> inline
    const TSFRBTreeNode<Key, Value>* TSFRBTree<Key, Value>::maximum() const
    {
      if (root_==0) return 0;
      else return root_->maximum();
    }




  template<class Key, class Value> inline
    const TSFRBTreeNode<Key, Value>* TSFRBTree<Key, Value>::pointerTo(const Key& key) const
    {
      TSFRBTreeNode<Key, Value>* tryPtr = root_;

      /*
        This code is generating a possibly uninitialized variable error at the
        highest optimization levels.
      */
      while (tryPtr != 0)
        {
          if (tryPtr->key_ == key) return tryPtr;
          if (key < tryPtr->key_)
            {
              tryPtr = tryPtr->left_;
            }
          else
            {
              tryPtr = tryPtr->right_;
            }
        }
      return tryPtr;
    }

  template<class Key, class Value> inline
    bool TSFRBTree<Key, Value>::containsKey(const Key& key) const
    {
      const TSFRBTreeNode<Key, Value>* ptr = pointerTo(key);
      if (ptr==0) return false;
      (Key&) cacheKey_ = key;
      (Value&) cacheValue_ = ptr->value_;
      return true;
    }

  template<class Key, class Value> inline
    const Value& TSFRBTree<Key, Value>::get(const Key& key) const
    {
      if (!containsKey(key)) TSFError::raise("TSFRBTree::get() key not found");
      return cacheValue_;
    }

  template<class Key, class Value> inline
    bool TSFRBTree<Key, Value>::get(const Key& key, Value& value) const
    {
      const TSFRBTreeNode<Key, Value>* ptr = pointerTo(key);
      if (ptr==0) return false;
      value = ptr->value_;
      return true;
    }

  template<class Key, class Value> inline
    ostream& TSFRBTree<Key, Value>::print(ostream& os, bool debug) const
    {
      if (root_!=0) root_->inorderPrint(os, debug);
      return os;
    }


  template<class Key, class Value> inline
    ostream& operator<<(ostream& os, const TSFRBTree<Key, Value>& t)
    {
      return t.print(os, false);
    }

  template<class Key, class Value> inline
    void TSFRBTree<Key, Value>::arrayify(TSFArray<Key>& keys,
                                         TSFArray<Value>& values) const
    {
      keys.resize(count_);
      values.resize(count_);

      const TSFRBTreeNode<Key, Value>* ptr = minimum();

      int i=0;
      while (ptr)
        {
          keys[i] = ptr->key_;
          values[i] = ptr->value_;
          ptr = ptr->successor();
          i++;
        }
      if (i != count_) TSFError::raise("TSFRBTree<Key, Value>::arrayify mismatch between"
                                       " count_ and number of elements");
    }

  template<class Key, class Value> inline
    TSFArray<Key> TSFRBTree<Key, Value>::arrayifyKeys() const
    {
      TSFArray<Key> keys(count_);

      const TSFRBTreeNode<Key, Value>* ptr = minimum();

      int i=0;
      while (ptr)
        {
          keys[i] = ptr->key_;
          ptr = ptr->successor();
          i++;
        }
      return keys;
    }

  template<class Key, class Value> inline
    TSFArray<Value> TSFRBTree<Key, Value>::arrayifyValues() const
    {
      TSFArray<Value> values(count_);

      const TSFRBTreeNode<Key, Value>* ptr = minimum();

      int i=0;
      while (ptr)
        {
          values[i] = ptr->value_;
          ptr = ptr->successor();
          i++;
        }
      return values;
    }





  //-------------------------------------------------------------------------
  //
  //                      Rotations


  template<class Key, class Value> inline
    void TSFRBTree<Key, Value>::rotateLeft(TSFRBTreeNode<Key, Value>* x)
    {
      // check for bad input; failure indicates a logic error in the
      // red-black balancing code.
      if (x==0) TSFError::raise("TSFRBTree::rotateLeft input error: x = nil. This should never happen, and indicates a bug in the tree balancing code.");
      if (x->right_ == 0) TSFError::raise("TSFRBTree::rotateLeft input error:  x->right = nil. This should never happen, and indicates a bug in the tree balancing code.");


      TSFRBTreeNode<Key, Value>* y = x->right_;

      // attach left subtree of y to right side of x.
      x->right_ = y->left_;
      // notify the former left subtree of y that it is now a subtree of x.
      if (y->left_ != 0) y->left_->parent_ = x;

      // link x's former parent to y.
      y->parent_ = x->parent_;

      // if x was root, set root to y.
      if (x->parent_ == 0)
        {
          root_ = y;
        }
      // otherwise, tell x's former parent about y.
      else if (x == x->parent_->left_)
        {
          x->parent_->left_ = y;
        }
      else
        {
          x->parent_->right_ = y;
        }

      // attach x to y's left
      y->left_ = x;
      // inform x that y is now its parent.
      x->parent_ = y;
    }




  template<class Key, class Value> inline
    void TSFRBTree<Key, Value>::rotateRight(TSFRBTreeNode<Key, Value>* x)
    {
      // check for bad input; failure indicates a logic error in the
      // red-black balancing code.
      if (x==0) TSFError::raise("TSFRBTree::rotateRight input error: x = nil. This should never happen, and indicates a bug in the tree balancing code.");
      if (x->left_ == 0) TSFError::raise("TSFRBTree::rotateRight input error:  x->left = nil. This should never happen, and indicates a bug in the tree balancing code.");


      TSFRBTreeNode<Key, Value>* y = x->left_;

      // attach right subtree of y to left side of x.
      x->left_ = y->right_;
      // notify the former right subtree of y that it is now a subtree of x.
      if (y->right_ != 0) y->right_->parent_ = x;

      // link x's former parent to y.
      y->parent_ = x->parent_;

      // if x was root, set root to y.
      if (x->parent_ == 0)
        {
          root_ = y;
        }
      // otherwise, tell x's former parent about y.
      else if (x == x->parent_->right_)
        {
          x->parent_->right_ = y;
        }
      else
        {
          x->parent_->left_ = y;
        }

      // attach x to y's right
      y->right_ = x;
      // inform x that y is now its parent.
      x->parent_ = y;
    }


}

#endif








