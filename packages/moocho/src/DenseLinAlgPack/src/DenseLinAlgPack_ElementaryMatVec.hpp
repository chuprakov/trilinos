// //////////////////////////////////////////////////////////////////////
// ElementaryMatVec.h

#ifndef ELEMENTARY_MAT_VEC_H
#define ELEMENTARY_MAT_VEC_H

#include "GenMatrixClass.h"

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
