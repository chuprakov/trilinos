#ifndef TSFSMARTPTR_H
#define TSFSMARTPTR_H

#include "TSFConfig.h"
#include <string>

// Solaris and IRIX CC can't handle inlines in some cases
#if defined (IRIX)
#define cond_inline 
#elif defined (SOLARIS)
#define cond_inline 
#else
#define cond_inline inline
#endif


namespace TSF
{
	using std::string;

	/**
	 * \ingroup Utilities
	 * Templated reference-counting pointer class. 
	 */

	template <class T>
		class TSFSmartPtr {
		public:
		cond_inline TSFSmartPtr(T* ptr = 0);
		cond_inline TSFSmartPtr(const TSFSmartPtr<T>& other);
		cond_inline TSFSmartPtr(T* ptr, int* refCount);
		~TSFSmartPtr();
		const TSFSmartPtr<T>& operator =(const TSFSmartPtr<T>& rhs);
		cond_inline T* operator ->();
		cond_inline bool isNull() const;
		cond_inline const T* operator ->() const;
		/** Read/write dereferencing */
		cond_inline T& operator *();
		/** Read-only dereferencing */
		cond_inline const T& operator *() const;
		cond_inline operator const T* () const;
		cond_inline bool isNonUnique() const;

		cond_inline int refCount() const;

		cond_inline void error(const string& msg) const ;

		protected:
		T* ptr_;
		int* refCount_;
	};

	void smartPtrError(const string& msg);

	template <class T> cond_inline
		TSFSmartPtr<T>::TSFSmartPtr(T* ptr)
		: ptr_(ptr),
		refCount_(0)
		{
			if (ptr_) {
				refCount_ = new int;
				if (refCount_ == 0)
					error("TSFSmartPtr::TSFSmartPtr(T* ptr) out of memory");
				*refCount_ = 1;
			}
		}

	template <class T> cond_inline
		TSFSmartPtr<T>::TSFSmartPtr(T* ptr, int* refCount)
		: ptr_(ptr),
		refCount_(refCount)
		{
			*refCount++;
		}

	template <class T> cond_inline
		TSFSmartPtr<T>::TSFSmartPtr(const TSFSmartPtr<T>& other)
		: ptr_(other.ptr_),
		refCount_(other.refCount_)
		{
			if (refCount_ != 0)
				++(*refCount_);
		}

	template <class T> 
		TSFSmartPtr<T>::~TSFSmartPtr() 
		{
			if (refCount_ != 0 && --(*refCount_) == 0) 
				{
					delete ptr_;
					ptr_ = 0;
					delete refCount_;
					refCount_ = 0;
				}
		}

	template <class T> 
		const TSFSmartPtr<T>& TSFSmartPtr<T>::operator =(const TSFSmartPtr<T>& rhs) 
	{
		if (ptr_ != rhs.ptr_) 
			{
				if (refCount_ != 0 && --(*refCount_) == 0) 
					{
						delete ptr_;
						delete refCount_;
					}
				ptr_ = rhs.ptr_;
				refCount_ = rhs.refCount_;
				if (refCount_ != 0)
					++(*refCount_);
			}
		return *this;
	}

	template <class T> cond_inline
		bool TSFSmartPtr<T>::isNull() const {
		return (ptr_ == 0);
	}

	template <class T> cond_inline
		T* TSFSmartPtr<T>::operator ->() {

		if (ptr_ == 0)
			error("TSFSmartPtr<T>::operator ->() on null pointer");
		return ptr_;
	}

	template <class T> cond_inline
		const T* TSFSmartPtr<T>::operator ->() const {
		if (ptr_ == 0)
			error("TSFSmartPtr<T>::operator ->() on null pointer");
		return ptr_;
	}

	template <class T> cond_inline
		T& TSFSmartPtr<T>::operator *() {
		if (ptr_ == 0)
			error("TSFSmartPtr<T>::operator *() on null pointer");
		return *ptr_;
	}

	template <class T> cond_inline
		const T& TSFSmartPtr<T>::operator *() const {
		if (ptr_ == 0)
			error("TSFSmartPtr<T>::operator *() on null pointer");
		return *ptr_;
	}

	template <class T> cond_inline
		TSFSmartPtr<T>::operator const T* () const {
		return ptr_;
	}

	template <class T> cond_inline
		bool TSFSmartPtr<T>::isNonUnique() const {
		return refCount_ == 0 ? false : *refCount_ != 1;
	}

	template <class T> cond_inline
		int TSFSmartPtr<T>::refCount() const {
		return refCount_ == 0 ? 0 : *refCount_;
	}
	
	template <class T>
		void TSFSmartPtr<T>::error(const string& msg) const 
		{
			smartPtrError(msg);
		}

	template <class T>
		string toString(const TSFSmartPtr<T>& sp) 
		{
			return "SmartPtr<" + TSF::toString(*sp) + ">";
		}
}

#undef cond_inline
#endif
