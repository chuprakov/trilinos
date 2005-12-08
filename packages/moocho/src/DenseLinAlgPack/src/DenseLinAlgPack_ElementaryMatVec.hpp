// //////////////////////////////////////////////////////////////////////
// DenseLinAlgPack_ElementaryMatVec.hpp
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

#include "DenseLinAlgPack_DMatrixClass.hpp"

namespace DenseLinAlgPack {

/// creates an n sized vector with all zeros accepts the ith element is one.
inline DVector e_vec(size_type i, size_type n) {
	DenseLinAlgPack::DVector v(0.0,n);
	v(i) = 1.0;
	return v;
}

/// creates an n x n identity matrix
inline DMatrix eye(size_type n) {
	DMatrix mat(0.0,n,n);
	mat.diag() = 1.0;
	return mat;
}

}	// namespace DenseLinAlgPack

#endif // ELEMENTARY_MAT_VEC_H
