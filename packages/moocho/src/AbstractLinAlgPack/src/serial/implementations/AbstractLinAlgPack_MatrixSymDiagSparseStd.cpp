// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymDiagonalSparseStd.cpp
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

#include "SparseLinAlgPack/include/MatrixSymDiagonalSparseStd.h"
#include "AbstractLinAlgPack/include/SpVectorOut.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"
#include "ThrowException.h"

namespace SparseLinAlgPack {

MatrixSymDiagonalSparseStd::MatrixSymDiagonalSparseStd( const SpVectorSlice& diag )
	: diag_(diag)
{}

void MatrixSymDiagonalSparseStd::initialize( const SpVectorSlice& diag )
{
	diag_ = diag;
}

// Overridden from MatrixWithOp

MatrixWithOp& MatrixSymDiagonalSparseStd::operator=(const MatrixWithOp& m)
{
	if(&m == this) return *this;	// assignment to self
	const MatrixSymDiagonalSparseStd
		*p_m = dynamic_cast<const MatrixSymDiagonalSparseStd*>(&m);
	if(p_m) {
		diag_ = p_m->diag_;
	}
	else {
		THROW_EXCEPTION(
			true, std::invalid_argument
			,"MatrixSymDiagonalSparseStd::operator=(const MatrixWithOp& m) : Error!"
			"The concrete type of m = \'" << typeid(m).name() << "\' is not a subclass of "
			"MatrixSymDiagonalSparseStd as expected"
			);
	}
	return *this;
}

// Overridden from MatrixDiagonalSparse

const SpVectorSlice MatrixSymDiagonalSparseStd::diag() const
{
	return diag_();
}


}	// end namespace SparseLinAlgPack
