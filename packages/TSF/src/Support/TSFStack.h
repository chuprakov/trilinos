#ifndef TSFSTACK_H
#define TSFSTACK_H

#include "TSFArray.h"

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
			TSFArray<T> list_;
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
			if (list_.length()==0) TSFError::raise("TSFStack<T>::get() on empty TSFStack");
            T rtn = list_[list_.length()-1];
			list_.remove(list_.length()-1);
			return rtn;
		}

	template<class T> inline T TSFStack<T>::peek() 
		{
			if (list_.length()==0) TSFError::raise("TSFStack<T>::get() on empty TSFStack");
			return list_[list_.length()-1];
		}


	template<class T> inline TSFArray<T> TSFStack<T>::arrayify()
		{
          return list_;
		}

	template<class T> inline string toString(const TSFStack<T>& stack)
		{
			return toString(list_);
		}



}
#endif





