// //////////////////////////////////////////////////////////////////////
// MatrixWithOpNonsingular.cpp
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

#include "AbstractLinAlgPack/include/MatrixWithOpNonsingular.h"

namespace AbstractLinAlgPack {

MatrixWithOpNonsingular::mat_mwons_mut_ptr_t
MatrixWithOpNonsingular::clone_mwons()
{
	return ReferenceCountingPack::null;
}

MatrixWithOpNonsingular::mat_mwons_ptr_t
MatrixWithOpNonsingular::clone_mwons() const
{
	return ReferenceCountingPack::null;
}

// Overridden from MatrixWithOp

MatrixWithOpNonsingular::mat_mut_ptr_t
MatrixWithOpNonsingular::clone()
{
	return clone_mwons();
}

MatrixWithOpNonsingular::mat_ptr_t
MatrixWithOpNonsingular::clone() const
{
	return clone_mwons();
}

// Overridden from MatrixNonsingular

MatrixWithOpNonsingular::mat_mns_mut_ptr_t
MatrixWithOpNonsingular::clone_mns()
{
	return clone_mwons();
}

MatrixWithOpNonsingular::mat_mns_ptr_t
MatrixWithOpNonsingular::clone_mns() const
{
	return clone_mwons();
}

}	// end namespace AbstractLinAlgPack
