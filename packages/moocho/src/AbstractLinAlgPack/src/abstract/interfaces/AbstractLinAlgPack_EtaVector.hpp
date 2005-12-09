// ///////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_EtaVector.hpp
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

#ifndef ETA_VECTOR_H
#define ETA_VECTOR_H

#include "AbstractLinAlgPack_SpVectorClass.hpp"

namespace AbstractLinAlgPack {

///
/** Create an eta vector (scaled by alpha = default 1).
  *
  * The created vector is of size n and has the single nonzero
  * element of eta(i) = alpha.
  * 
  * The default constructor and assignment functions are not
  * allowed.
  */
class EtaVector {
public:
	
	typedef SpVectorSlice::element_type		ele_t;


	///
	EtaVector( ele_t::index_type i, size_type n, ele_t::value_type alpha = 1.0 )
		: ele_(i,alpha), sv_(&ele_,1,0,n,true)
	{}

	/// Implicit conversion to a SpVectorSlice object.
	operator const SpVectorSlice() const
	{
		return sv_;
	}

	/// Explicit conversion to a SpVectorSlice object.
	const SpVectorSlice& operator()() const
	{
		return sv_;
	}

private:
	ele_t			ele_;
	SpVectorSlice	sv_;

	// not defined and not to be called
	EtaVector();
	EtaVector& operator=(const EtaVector&);

};	// end class EtaVector


}	// namespace AbstractLinAlgPack

#endif	// ETA_VECTOR_H
