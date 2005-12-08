// //////////////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_IVector.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef IVECTOR_H
#define IVECTOR_H

#include <assert.h>

#include <valarray>

#include "DenseLinAlgPack_Types.hpp"
#include "Teuchos_TestForException.hpp"

namespace DenseLinAlgPack {
///
/* * Fortran compatable integer vector for holding the  pivot information for
 * the elements of a vector, or the rows or columns of a matrix.
 */
class IVector : public std::valarray<DenseLinAlgPack::size_type> {
public:

	// STL typedefs
	typedef DenseLinAlgPack::index_type		value_type;
	typedef DenseLinAlgPack::size_type		size_type;
	typedef value_type&					reference;
	typedef const value_type&			const_reference;
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

	/// Resize on assignment
	IVector& operator=(const IVector&);

	/// 1-based element access (range checked if _DEBUG is defined)
	reference operator()(size_type i);
	/// 1-based element access (range checked if _DEBUG is defined)
	const_reference operator()(size_type i) const;

	/// STL iterator
	iterator begin();
	/// STL iterator
	const_iterator begin() const;
	/// STL iterator
	iterator end();
	/// STL iterator
	const_iterator end() const;

}; // end class IVector

// Inline definitions

inline IVector::IVector() : std::valarray<size_type>()
{}

inline IVector::IVector(size_type n) : std::valarray<size_type>(n)
{}

inline IVector::IVector(const value_type& val, size_type n) : std::valarray<size_type>(val,n)
{}

inline IVector::IVector(const value_type* p, size_type n) : std::valarray<size_type>(p,n)
{}

inline IVector& IVector::operator=(const IVector& iv)
{
	this->resize(iv.size());
	std::valarray<DenseLinAlgPack::size_type>::operator=(iv);
	return *this;
}

inline IVector::reference IVector::operator()(size_type i)
{
#ifdef _DEBUG
	assert( 1 <= i && i <= static_cast<size_type>(size()) );
#endif
	return operator[](i-1);
}

inline IVector::const_reference IVector::operator()(size_type i) const
{
#ifdef _DEBUG
	assert( 1 <= i && i <= static_cast<size_type>(size()) );
#endif
	return const_cast<IVector*>(this)->operator[](i-1);
}

inline IVector::iterator IVector::begin()
{	return &operator[](0); }

inline IVector::const_iterator IVector::begin() const
{	return &(const_cast<IVector*>(this)->operator[](0)); }

inline IVector::iterator IVector::end()
{	return begin() + size(); }

inline IVector::const_iterator IVector::end() const
{	return begin() + size(); }

}	// end namespace DenseLinAlgPack

#endif // IVECTOR_H
