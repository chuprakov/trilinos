// /////////////////////////////////////////////////////
// ReleaseResource_ref_count_ptr.hpp
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

#ifndef RELEASE_RESOURCE_REF_COUNT_PTR_H
#define RELEASE_RESOURCE_REF_COUNT_PTR_H

#include "ReleaseResource.hpp"
#include "ref_count_ptr.hpp"

namespace MemMngPack {

///
/** Template class that implements ReleaseResource interface for
 * a ref_count_ptr<T> object.
 *
 * Note that ~ReleaseResource_ref_count_ptr() does not need to be
 * implemented since the compiler generated version will already
 * be correct.
 */
template <class T>
class ReleaseResource_ref_count_ptr : public ReleaseResource {
public:

	///
	typedef MemMngPack::ref_count_ptr<T>   ptr_t;

	/// Just give public access to pointer
	ptr_t  ptr;

	/// Construct from a pointer
	ReleaseResource_ref_count_ptr(const ptr_t& ptr);

	// ////////////////////////////////////
	// Overriddend from ReleaseResource

	///
	bool resource_is_bound() const;

}; // end class ReleaseResource_ref_count_ptr

// ////////////////////////////////
// Inline function definitions

template <class T>
inline
ReleaseResource_ref_count_ptr<T>::ReleaseResource_ref_count_ptr(const ptr_t& p)
	: ptr(p)
{}

// ///////////////////////////////
// Template function definitions

template <class T>
bool ReleaseResource_ref_count_ptr<T>::resource_is_bound() const
{
	return ptr.get() != 0;
}

} // end namespace MemMngPack

#endif // RELEASE_RESOURCE_REF_COUNT_PTR_H
