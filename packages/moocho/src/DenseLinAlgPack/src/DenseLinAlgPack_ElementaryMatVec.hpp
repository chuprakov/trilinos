// //////////////////////////////////////////////////////////////////////
// ElementaryMatVec.hpp
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

#ifndef ELEMENTARY_MAT_VEC_H
#define ELEMENTARY_MAT_VEC_H

#include "GenMatrixClass.hpp"

namespace LinAlgPack {

/// creates an n sized vector with all zeros accepts the ith element is one.
inline Vector e_vec(size_type i, size_type n) {
	LinAlgPack::Vector v(0.0,n);
	v(i) = 1.0;
	return v;
}

/// creates an n x n identity matrix
inline GenMatrix eye(size_type n) {
	GenMatrix mat(0.0,n,n);
	mat.diag() = 1.0;
	return mat;
}

}	// namespace LinAlgPack

#endif // ELEMENTARY_MAT_VEC_H
