// //////////////////////////////////////////////////////////////////////////////////
// MatrixSymDiagSparseStd.cpp
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

#include "AbstractLinAlgPack/src/serial/implementations/MatrixSymDiagSparseStd.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/SpVectorOut.hpp"
#include "DenseLinAlgPack/src/DVectorClass.hpp"
#include "DenseLinAlgPack/src/DenseLinAlgPackAssertOp.hpp"
#include "ThrowException.hpp"

namespace AbstractLinAlgPack {

MatrixSymDiagSparseStd::MatrixSymDiagSparseStd( const SpVectorSlice& diag )
	: diag_(diag)
{}

void MatrixSymDiagSparseStd::initialize( const SpVectorSlice& diag )
{
	diag_ = diag;
}

// Overridden from MatrixOp

MatrixOp& MatrixSymDiagSparseStd::operator=(const MatrixOp& m)
{
	if(&m == this) return *this;	// assignment to self
	const MatrixSymDiagSparseStd
		*p_m = dynamic_cast<const MatrixSymDiagSparseStd*>(&m);
	if(p_m) {
		diag_ = p_m->diag_;
	}
	else {
		THROW_EXCEPTION(
			true, std::invalid_argument
			,"MatrixSymDiagSparseStd::operator=(const MatrixOp& m) : Error!"
			"The concrete type of m = \'" << typeid(m).name() << "\' is not a subclass of "
			"MatrixSymDiagSparseStd as expected"
			);
	}
	return *this;
}

// Overridden from MatrixDiagonalSparse

const SpVectorSlice MatrixSymDiagSparseStd::diag() const
{
	return diag_();
}


}	// end namespace AbstractLinAlgPack
