#ifndef TSFSTACK_H
#define TSFSTACK_H

#include "TSFLinkedList.h"

namespace TSF
{
	using std::string;

	/**
	 * \ingroup Containers
	 * Templated LIFO stack.
	 * @author Kevin Long
	 */


	template<class T>
		class TSFStack
		{
		public:
			/** Construct an empty stack */
			inline TSFStack();
	
			/** push an element onto the top of the stack */
			inline void push(const T& data);

			/** pop the top element from the top of the stack */
			inline T pop();

			/** peek at the top element */
			inline T peek();

			/** get the number of elements in the stack */
			inline int size() const {return list_.length();}

			/** read elements into an Array */
			inline TSFArray<T> arrayify();

		private:
			TSFLinkedList<T> list_;
		};


	// create an empty list
	template<class T> inline TSFStack<T>::TSFStack()
		: list_()
		{;}

	// put a new entry at the beginning of the list

	template<class T> inline void TSFStack<T>::push(const T& data)
		{
			list_.append(data);
		}

	// return last entry from list, then unhook it from the list.

	template<class T> inline T TSFStack<T>::pop()
		{
			T rtn;
			if (!list_.getLastNode(rtn)) TSFError::raise("TSFStack<T>::get() on empty TSFStack");
			list_.removeLast();
			return rtn;
		}

	template<class T> inline T TSFStack<T>::peek() 
		{
			T rtn;
			if (!list_.getLastNode(rtn)) TSFError::raise("TSFStack<T>::get() on empty TSFStack");
			return rtn;
		}


	template<class T> inline TSFArray<T> TSFStack<T>::arrayify()
		{
			return ::arrayify(list_);
		}

	template<class T> inline string toString(const TSFStack<T>& stack)
		{
			return toString(stack.arrayify());
		}



}
#endif





