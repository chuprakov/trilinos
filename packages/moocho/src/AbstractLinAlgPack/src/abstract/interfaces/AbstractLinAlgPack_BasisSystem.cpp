// //////////////////////////////////////////////////////////////
// BasisSystem.cpp
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

#include "AbstractLinAlgPack/include/BasisSystem.h"
#include "Range1D.h"

namespace AbstractLinAlgPack {

Range1D BasisSystem::equ_decomp() const
{
	return Range1D(1,this->var_dep().size());
}

Range1D BasisSystem::equ_undecomp() const
{
	return Range1D::Invalid;
}

Range1D BasisSystem::inequ_decomp() const
{
	return Range1D::Invalid;
}

Range1D BasisSystem::inequ_undecomp() const
{
	return Range1D::Invalid;
}

} // end namespace AbstractLinAlgPack
