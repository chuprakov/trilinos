#ifndef TSFQUEUE_H
#define TSFQUEUE_H

#include "TSFConfig.h"
#include "TSFLinkedList.h"

namespace TSF
{

	/**
 * \ingroup Containers
 * Circular linked-list implementation of a FIFO queue.
 * 
 * @author Kevin Long
*/


	template<class T>
		class TSFQueue
		{
		public:
			/** Create an empty queue */
			TSFQueue();
	
			/** Put an entry at the back of the queue */
			void put(const T& data);
			/** Get the entry from the front of the queue */
			T get();

			/** determine if the queue is empty */
			bool isEmpty() const ;

			/** Read elements into an TSFArray */
			TSFArray<T> arrayify();

			string toString();
		private:
			TSFLinkedList<T> list_;
		};

	// create an empty list
	template<class T> TSFQueue<T>::TSFQueue()
		: list_()
		{;}

	// put a new entry at the beginning of the list

	template<class T> void TSFQueue<T>::put(const T& data)
		{
			list_.prepend(data);
		}

	// return last entry from list, then unhook it from the list.

	template<class T> T TSFQueue<T>::get()
		{
			T rtn;
			if (!list_.getLastNode(rtn)) TSFError::raise("TSFQueue<T>::get() on empty queue");
			list_.removeLast();
			return rtn;
		}

	template<class T> bool TSFQueue<T>::isEmpty() const
		{
			return list_.length()==0;
		}

	template<class T> inline string TSFQueue<T>::toString()
		{
			TSFLinkedList<T> tmp = list_;
			TSFArray<T> a = ::arrayify(tmp);
			return toString(a);
		}

	template<class T> inline string toString(const TSFQueue<T>& q)
		{
			return q.toString();
		}

	template<class T> inline TSFArray<T> TSFQueue<T>::arrayify()
		{
			return ::arrayify(list_);
		}

	template<class T> inline ostream& operator<<(ostream& os, const TSFQueue<T>& q)
		{
			return os << toString(q);
		}


}

#endif





