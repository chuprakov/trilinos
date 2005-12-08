// /////////////////////////////////////////////////////////////////////////////
// PermIn.cpp
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

#include "DenseLinAlgPack_PermIn.hpp"
#include "DenseLinAlgPack_IVector.hpp"

std::istream& DenseLinAlgPack::operator>>(std::istream& istrm, IVector& perm) {
	size_type size;
	istrm >> size;
	perm.resize(size);

	IVector::iterator		itr_perm		= perm.begin(),
							itr_perm_end	= perm.end();
	for(;itr_perm != itr_perm_end;)
		istrm >> *itr_perm++;
	return istrm;
}
