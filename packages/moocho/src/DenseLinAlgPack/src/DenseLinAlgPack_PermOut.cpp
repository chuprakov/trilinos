// /////////////////////////////////////////////////////////////////////////////
// PermOut.cpp
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

#include <iomanip>

#include "DenseLinAlgPack_PermOut.hpp"
#include "DenseLinAlgPack_IVector.hpp"

std::ostream& DenseLinAlgPack::operator<<(std::ostream& o, const IVector& perm) {
	int w = o.width(0) - 1; // get the set width
	o << perm.size() << "\n";
	IVector::const_iterator
		itr_perm		= perm.begin(),
		itr_perm_end	= perm.end();
	for(;itr_perm != itr_perm_end;)
		o << std::setw(w) << *itr_perm++ << ' ';
	o << std::endl;
	return o;
}
