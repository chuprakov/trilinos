// //////////////////////////////////////////////////////////////////////////////
// IVector.h
//
// Integer vector used for holding pivot information
//

#ifndef IVECTOR_H
#define IVECTOR_H

#include <valarray>

#include "LinAlgPackTypes.h"

namespace LinAlgPack {
///
/** Fortran compatable integer vector for holding the 
  * pivot information for the elements of a vector, or the rows or columns
  * of a matrix.
  */
class IVector : public std::valarray<LinAlgPack::size_type> {
public:
	// STL typedefs
	typedef LinAlgPack::size_type		value_type;
	typedef LinAlgPack::size_type		size_type;
	typedef value_type&					reference;
	typedef value_type*					iterator;
	typedef const value_type*			const_iterator;
	typedef std::valarray<size_type>	valarray;

	// constructors

	///
	IVector();
	///
	IVector(size_type n);
	///
	IVector(const value_type& val, size_type n);
	///
	IVector(const value_type* p, size_type n);

	/// 1-based element access
	reference operator()(size_type i);
	/// 1-based element access
	value_type operator()(size_type i) const;

	/// STL iterator
	iterator begin();
	/// STL iterator
	const_iterator begin() const;
	/// STL iterator
	iterator end();
	/// STL iterator
	const_iterator end() const;
};

// Inline definitions

inline IVector::IVector() : std::valarray<size_type>()
{}

inline IVector::IVector(size_type n) : std::valarray<size_type>(n)
{}

inline IVector::IVector(const value_type& val, size_type n) : std::valarray<size_type>(val,n)
{}

inline IVector::IVector(const value_type* p, size_type n) : std::valarray<size_type>(p,n)
{}

inline IVector::reference IVector::operator()(size_type i)
{	return operator[](i-1);	}

inline IVector::value_type IVector::operator()(size_type i) const
{	return operator[](i-1);		}

inline IVector::iterator IVector::begin()
{	return &operator[](0); }

inline IVector::const_iterator IVector::begin() const
{	return &(const_cast<IVector*>(this)->operator[](0)); }

inline IVector::iterator IVector::end()
{	return begin() + size(); }

inline IVector::const_iterator IVector::end() const
{	return begin() + size(); }

}	// end namespace LinAlgPack

#endif