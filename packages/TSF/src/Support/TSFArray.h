#ifndef TSFARRAY_H
#define TSFARRAY_H

#include <iostream>
#include <string>
#include "TSFError.h"
#include "TSFUtils.h"

using namespace std;
namespace TSF
{

	/** \ingroup Containers 
	 * TSFArray is a templated container class similar to the STL vector, but with
	 * a more lightweight implementation and API. 
	 */ 
	template<class T>
		class TSFArray 
		{
		public:
			/** Empty ctor */
			TSFArray();

			/** Allocate an array with n elements */
			TSFArray(int n); 

			/** Allocate n elements, and fill with value */
			TSFArray(int n, const T& t); 

			/** dtor deletes the contents */
			~TSFArray();

			/** copy ctor makes a deep copy of all elements */
			TSFArray(const TSFArray<T>& other);  

			/** assignment operator makes a deep copy of all elements */
			const TSFArray<T>& operator=(const TSFArray<T>& other); 

			/** 
			 * return number of elements.
			 */
			int length() const ;
			int size() const ;

			/** change size */
			void resize(int newN);
			/** preallocate space */
			void reserve(int n); 
			/** get preallocated space */
			int  capacity() const ;

			/** Stick a new entry at the end of the array. Resize to allow
			 * space for the new entry 
			 */
			TSFArray<T>& append(const T& rhs);
	
			/** 
			 * Read/Write access to a the i-th element. 
			 */
			T& operator[](int i);
			/** 
			 * Read-only access to a the i-th element. 
			 */
			const T& operator[](int i) const;

			/** 
			 * Remove the i-th element. Subsequent elements are bumped back 
			 * to fill in the hole.
			 */
			void remove(int i);

			/** */
			string toString() const ; 
	
		private:
			/** C array of contents */
			T* data_;     
			/** number of active elements */
			int len_;     
			/** number of elements that have been allocated */
			int capacity_; 

			/** check for a bounds violation if HAVE_ARRAY_BOUNDSCHECK has been
			 * defined as 1. */
			void indexCheckCrash(int i) const;
		};

	/** \relates TSFArray */
	template<class T> ostream& operator<<(ostream& os, const TSFArray<T>& array); 

	/** \relates TSFArray */
	template<class T> int hashCode(const TSFArray<T>& array);

	/** \relates TSFArray */
	template<class T> string toString(const TSFArray<T>& array);

	template<class T> inline TSFArray<T>::TSFArray()
		: data_(0),
    len_(0), 
    capacity_(0)
		{
		}



	template<class T> inline TSFArray<T>::TSFArray(int n)
		: data_(0),
    len_(n), 
    capacity_(n)
		{
			if (len_ < 0)
				TSFError::raise("Negative length passed to TSFArray<T>::TSFArray(int n)");
			if (len_ > 0) 
				{
					data_ = new T [len_];
					if (!data_)
						TSFError::raise("TSFArray constructor out of memory");
				}
		}


	template<class T> inline TSFArray<T>::TSFArray(int n, const T& t)
		: data_(0),
    len_(n), 
    capacity_(n)
		{
			if (len_ < 0)
				TSFError::raise("Negative length passed to TSFArray<T>::TSFArray(int n)");
			if (len_ > 0) 
				{
					data_ = new T [len_];
					if (!data_)
						TSFError::raise("TSFArray constructor out of memory");
				}
			for (int i = 0; i < len_; ++ i) 
				data_[i] = t;
		}


	template<class T> inline  TSFArray<T>::TSFArray(const TSFArray<T>& arr)
		: data_(0),
    len_(arr.len_), 
    capacity_(arr.len_)
		{
			if (len_ > 0) 
				{
					data_ = new T [capacity_];
					if (!data_)
						TSFError::raise("TSFArray constructor out of memory");
					for (int i = 0; i < len_; ++i)
						data_[i] = arr.data_[i];
				}
		}

	template<class T> inline TSFArray<T>::~TSFArray() 
		{
			delete [] data_;
		}

	template<class T> inline 
		const TSFArray<T>& TSFArray<T>::operator=(const TSFArray<T>& arr) 
	{
		if (this != &arr) 
			{ // don't bother to assign if they're already identical
				if (capacity_ < arr.len_) 
					{ //If the reserved space is too small to hold arr
						delete [] data_;
						data_ = 0;
						capacity_ = arr.len_;
						if (capacity_ > 0) {
							data_ = new T[capacity_];
							if (!data_)
								TSFError::raise("TSFArray constructor out of memory");
						}
					}
				len_ = arr.len_;
				for (int i = 0; i < len_; ++i)
					data_[i] = arr[i];
			}
		return *this;
	}



	inline int bump(int start, int finish) 
		{
			if (finish < 10000)
				{
					if (start == 0)
						start = 1;
					while (start < finish)
						start *= 2;
					return start;
				}
			else
				{
					while (start < finish)
						{
							start += 10000;
						}
					return start;
				}
		}


	template<class T> inline
		void TSFArray<T>::resize(int newN) {
		if (len_ != newN) { // do not do anything if new size is not different
			if (newN < 0)
				TSFError::raise("Negative length passed to TSFArray<T>::resize(int newN)");
			if(newN > capacity_)
				reserve(bump(capacity_, newN));
			len_ = newN;
		}
	}




	template<class T> inline
		void TSFArray<T>::reserve(int N){
		if(capacity_ != N){
			if(N < 0){
				TSFError::raise("Negative length passed to TSFArray<T>::reserve(int N)");
			}
			if(N < len_){ len_ = N;}
			capacity_ = N;
			T* oldData = data_;
			data_ = 0;
			data_ = new T [capacity_];
			if (!data_)
				TSFError::raise("TSFArray<T>::reserve(int N) out of memory");
			for (int i = 0; i < len_; i++)
				data_[i] = oldData[i];
			delete [] oldData;
		}
	}

	template<class T> inline
		int TSFArray<T>::capacity() const{
		return capacity_;
	}



	template<class T> inline
		TSFArray<T>& TSFArray<T>::append(const T& rhs) 
		{
			resize(len_+1);
			data_[len_-1] = rhs;
			return *this;
		}



	template<class T> inline
		T& TSFArray<T>::operator[](int i) {
#if HAVE_ARRAY_BOUNDSCHECK
		indexCheckCrash(i);
#endif
		return data_[i];
	}

	template<class T> inline
		const T& TSFArray<T>::operator[](int i) const {
#if HAVE_ARRAY_BOUNDSCHECK
		indexCheckCrash(i);
#endif
		return data_[i];
	}

	template<class T> inline
		void TSFArray<T>::remove(int i)
		{
#if HAVE_ARRAY_BOUNDSCHECK
			indexCheckCrash(i);
#endif
			for (int j=i+1; j<length(); j++)
				{
					data_[j-1] = data_[j];
				}
			data_[len_-1] = T();
			len_--;
		}

	template<class T> inline
		int TSFArray<T>::length() const {
		return len_;
	}

	template<class T> inline
		int TSFArray<T>::size() const {
		return len_;
	}


	template<class T> inline
		void TSFArray<T>::indexCheckCrash(int i) const 
		{
			if (i<0 || i>=len_)
				TSFError::raise("TSFArray<T> index=" + TSF::toString(i) + " out of range [0,"
												+ TSF::toString(len_) + ")");
		}

	template<class T>
		TSFArray<T> sliceTSFArray(const TSFArray<TSFArray<T> >& array, int index)
		{
			TSFArray<T> rtn(array.length());

			for (int i=0; i<array.length(); i++)
				{
					rtn[i] = array[i][index];
				}
			return rtn;
		}

	
	// print in form (), (1), or (1,2)
	template<class T> inline ostream& operator<<(ostream& os, const TSFArray<T>& array) 
		{
			return os << toString(array);
		}

	template<class T> inline int hashCode(const TSFArray<T>& array)
		{
			int rtn = hashCode(len_);
			for (int i=0; i<a.length(); i++)
				{
					rtn += hashCode(a[i]);
				}
			return rtn;
		}

	template<class T> inline string TSFArray<T>::toString() const 
		{
			string rtn = "{";

			for (int i=0; i<length(); i++)
				{
					rtn += TSF::toString(data_[i]);
					if (i<length()-1) rtn += ", ";
				}
			rtn += "}";

			return rtn;
		}

	template<class T> inline string toString(const TSFArray<T>& array)
		{
			return array.toString();
		}

	// utilities for building small arrays

	template<class T> inline
		TSFArray<T> tuple(const T& a)
		{
			TSFArray<T> rtn(1, a);
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b)
		{
			TSFArray<T> rtn(2);
			rtn[0] = a;
			rtn[1] = b;
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c)
		{
			TSFArray<T> rtn(3);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c, const T& d)
		{
			TSFArray<T> rtn(4);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			rtn[3] = d;
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
		{
			TSFArray<T> rtn(5);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			rtn[3] = d;
			rtn[4] = e;
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
											const T& f)
		{
			TSFArray<T> rtn(6);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			rtn[3] = d;
			rtn[4] = e;
			rtn[5] = f;
			return rtn;
		}

	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
											const T& f, const T& g)
		{
			TSFArray<T> rtn(7);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			rtn[3] = d;
			rtn[4] = e;
			rtn[5] = f;
			rtn[6] = g;
			return rtn;
		}


	template<class T> inline
		TSFArray<T> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
											const T& f, const T& g, const T& h)
		{
			TSFArray<T> rtn(8);
			rtn[0] = a;
			rtn[1] = b;
			rtn[2] = c;
			rtn[3] = d;
			rtn[4] = e;
			rtn[5] = f;
			rtn[6] = g;
			rtn[7] = h;
			return rtn;
		}

}
	
#endif

