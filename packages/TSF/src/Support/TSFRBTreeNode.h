#ifndef TSFRBTREENODE_H
#define TSFRBTREENODE_H

#include "TSFConfig.h"
#include "TSFError.h"
#include <iostream>

namespace TSF
{
	using std::ostream;

	/** */
	enum TSFRBColor {TSFRBRed, TSFRBBlack};

	/** \ingroup Containers A single node in a red-black tree. Each node
		 contains a key of type Key, a data field of type Value, a color
		 TSFRBRed | TSFRBBlack, and pointers to its parent and left/right child
		 nodes.

	*/

	template<class Key, class Value> class TSFRBTreeNode
		{
		public:
			/** ctor */
			inline TSFRBTreeNode(const Key& key, const Value& value,
												TSFRBTreeNode<Key, Value>* parent, TSFRBColor color);
			inline ~TSFRBTreeNode();
			inline TSFRBTreeNode<Key, Value>* clone() const ;

			// printing, with variable verboseness.
			inline void inorderPrint(ostream& os, bool debug) const ;
			// verbose printing.
			inline void debugPrint(ostream& os) const ;
			// laconic printing.
			inline void print(ostream& os) const ;

			// traversal functions.
			inline const TSFRBTreeNode<Key, Value>* minimum() const ;
			inline const TSFRBTreeNode<Key, Value>* maximum() const ;
			inline const TSFRBTreeNode<Key, Value>* successor() const ;

			// some geneological operations that will be needed in 
			// balancing.
			inline bool isLeftChild() const ;
			inline TSFRBTreeNode<Key, Value>* grandparent() const ;


			// data members are public. They should be used only by TSFRBTRee.
			Key key_;
			Value value_;
			TSFRBColor color_;
			TSFRBTreeNode<Key, Value>* parent_;
			TSFRBTreeNode<Key, Value>* left_;
			TSFRBTreeNode<Key, Value>* right_;


		private:
			// copy ctor and assignment are private. 
			// The copy ctor should be called only by clone().
			// The assignment op should never be called. It is left undefined
			// and made private to prevent accidental use.
			inline TSFRBTreeNode(const TSFRBTreeNode<Key,Value>& other);
			inline const TSFRBTreeNode<Key,Value>& operator=
				(const TSFRBTreeNode<Key,Value>& other);
		};


	// Constructor: insert (key, value), assign a color, and link to a parent.
	template<class Key, class Value> inline
		TSFRBTreeNode<Key, Value>::TSFRBTreeNode(const Key& key, const Value& value,
																			 TSFRBTreeNode<Key, Value>* parent, 
																			 TSFRBColor color)
		: key_(key),
		value_(value),
		color_(color),
		parent_(parent),
		left_(0),
		right_(0)
		{
			;
		}

	// Make a copy by copying the left and right subtrees as well as the data
	// fields. Don't set the parent here; the calling code is responsible for
	// registering the new parent with this node.

	template<class Key, class Value> inline
		TSFRBTreeNode<Key, Value>::TSFRBTreeNode(const TSFRBTreeNode<Key,Value>& other)
		: key_(other.key_),
		
		value_(other.value_),
	  color_(other.color_),
	  parent_(0),
	  left_(0),
		right_(0)
		{
			if (other.left_ != 0)
				{
					left_ = other.left_->clone();
				}
			if (other.right_ != 0)
				{
					right_ = other.right_->clone();
				}
			// don't clone the parent. It will have to be set later.
		}


	// destructor: just delete the left and right subtrees.

	template<class Key, class Value> inline
		TSFRBTreeNode<Key, Value>::~TSFRBTreeNode() 
		{
			delete left_;
			delete right_;
		}


	// Make a copy of the node and its subtrees.
	// Register the copy node as the parent of the copy subtrees.

	template<class Key, class Value> inline
		TSFRBTreeNode<Key,Value>* TSFRBTreeNode<Key, Value>::clone() const 
		{
			TSFRBTreeNode<Key,Value>* result = new TSFRBTreeNode<Key,Value>(*this);
			if (result == 0)
				TSFError::raise("TSFRBTreeNode<Key, Value>::clone() out of memory");
			if (result->left_ != 0)
				{
					result->left_->parent_ = result;
				}
			if (result->right_ != 0)
				{
					result->right_->parent_ = result;
				}
	
			return result;
		}

	// TSFRBTreeNode* mininum()
	// Find the pointer to the leftmost child. This is easy: at every node, 
	// just walk left.

	template<class Key, class Value> inline 
		const TSFRBTreeNode<Key, Value>* TSFRBTreeNode<Key, Value>::minimum() const 
		{
			const TSFRBTreeNode<Key, Value>* minPtr = this;
			while (minPtr->left_ != 0)
				{
					minPtr = minPtr->left_;
				}
			return minPtr;
		}


	// TSFRBTreeNode* maximum()
	// Find the pointer to the rightmost child. This is easy: at every node, 
	// just walk right.

	template<class Key, class Value> inline 
		const TSFRBTreeNode<Key, Value>* TSFRBTreeNode<Key, Value>::maximum() const 
		{
			const TSFRBTreeNode<Key, Value>* maxPtr = this;
			while (maxPtr->right_ != 0)
				{
					maxPtr = maxPtr->right_;
				}
			return maxPtr;
		}


	// TSFRBTreeNode* successor()
	// Find the pointer to the next-highest-keyed node.
	// There are two possibilities:
	// (1) There is a right child tree. Then, find the minimum
	// of the right child subtree.
	// (2) There is no right child tree. In that case, walk back up the tree
	// until we find the node of which this branch is a left subtree.
	// (Refer to _Algorithms" figure 13.2 for a picture of this case).

	template<class Key, class Value> inline 
		const TSFRBTreeNode<Key, Value>* TSFRBTreeNode<Key, Value>::successor() const 
		{
			if (right_ != 0)
				{
					return right_->minimum();
				}
			const TSFRBTreeNode<Key, Value>* y = parent_;
			const TSFRBTreeNode<Key, Value>* x = this;
	
			while (y != 0 && x==y->right_)
				{
					x = y;
					y = y->parent_;
				}
			return y;
		}

	template<class Key, class Value> inline 
		bool TSFRBTreeNode<Key, Value>::isLeftChild() const
		{
			if (parent_->left_ == this) return true;
			return false;
		}

	template<class Key, class Value> inline 
		TSFRBTreeNode<Key, Value>* TSFRBTreeNode<Key, Value>::grandparent() const 
		{
			return parent_->parent_;
		}


	// Print nodes in order of keys. 
	// If debug==true, do a print of everything including pointers.
	// If debug==false, just print (key, value) pairs.

	template<class Key, class Value> inline
		void TSFRBTreeNode<Key, Value>::inorderPrint(ostream& os, bool debug) const 
		{
			if (left_) 
				{
					left_->inorderPrint(os, debug);
				}
			if (debug) debugPrint(os);
			else print(os);
	
			if (right_) 
				{
					if (!debug) os << ", ";
					right_->inorderPrint(os, debug);
				}
		}



	// Print everything including pointers.

	template<class Key, class Value> inline
		void TSFRBTreeNode<Key, Value>::debugPrint(ostream& os) const 
		{
			os << "P=" << parent_ 
				 << " This=" << this 
				 << " L=" << left_ 
				 << " R=" << right_
				 << " c=" << (int) color_ 
				 << " data=(" << key_ 
				 << ", " 
				 << value_ << ")" << endl;
		}

	// Print just the key and value

	template<class Key, class Value> inline
		void TSFRBTreeNode<Key, Value>::print(ostream& os) const 
		{
			os << "(" << key_ 
				 << ", " 
				 << value_ << ")";
		}


}

#endif
