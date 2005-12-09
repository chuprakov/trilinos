// /////////////////////////////////////////////////////////////////////////////
// DecompositionSystem.cpp
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

#include "ConstrainedOptPack_DecompositionSystem.hpp"

namespace ConstrainedOptPack {

size_type DecompositionSystem::n() const
{
	return this->space_range()->dim() + this->space_null()->dim();
}

size_type DecompositionSystem::r() const
{
	return this->space_range()->dim();
}

Range1D DecompositionSystem::equ_decomp() const
{
	const size_type
		r = this->r();
	return r ? Range1D(1,this->r()) : Range1D::Invalid;
}
	
Range1D DecompositionSystem::equ_undecomp() const
{
	const size_type
		r = this->r(),
		m = this->m();
	return m > r ? Range1D(r+1,m) : Range1D::Invalid;
}

}	// end namespace ConstrainedOptPack
