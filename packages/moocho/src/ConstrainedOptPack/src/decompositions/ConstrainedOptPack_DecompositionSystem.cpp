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

#include "ConstrainedOptimizationPack/include/DecompositionSystem.h"

namespace ConstrainedOptimizationPack {

Range1D DecompositionSystem::con_decomp() const
{
	return Range1D(1,this->r());
}
	
Range1D DecompositionSystem::con_undecomp() const
{
	const size_type
		r = this->r(),
		m = this->m();
	return m > r ? Range1D(r+1,m) : Range1D::Invalid;
}

}	// end namespace ConstrainedOptimizationPack
