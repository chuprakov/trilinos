// //////////////////////////////////////////////////////////////////////
// MatrixSymWithOpNonsingular.cpp
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

#include "AbstractLinAlgPack/src/MatrixSymWithOpNonsingular.hpp"

namespace AbstractLinAlgPack {

MatrixSymWithOpNonsingular::mat_mswons_mut_ptr_t
MatrixSymWithOpNonsingular::clone_mswons()
{
	return MemMngPack::null;
}

MatrixSymWithOpNonsingular::mat_mswons_ptr_t
MatrixSymWithOpNonsingular::clone_mswons() const
{
	return MemMngPack::null;
}

// Overridden from MatrixWithOp

MatrixSymWithOpNonsingular::mat_mut_ptr_t
MatrixSymWithOpNonsingular::clone()
{
	return clone_mswons();
}

MatrixSymWithOpNonsingular::mat_ptr_t
MatrixSymWithOpNonsingular::clone() const
{
	return clone_mswons();
}

// Overridden from MatrixNonsingular

MatrixSymWithOpNonsingular::mat_mns_mut_ptr_t
MatrixSymWithOpNonsingular::clone_mns()
{
	return clone_mswons();
}

MatrixSymWithOpNonsingular::mat_mns_ptr_t
MatrixSymWithOpNonsingular::clone_mns() const
{
	return clone_mswons();
}

// Overridden from MatrixSymWithOp

MatrixSymWithOpNonsingular::mat_mswo_mut_ptr_t
MatrixSymWithOpNonsingular::clone_mswo()
{
	return clone_mswons();
}

MatrixSymWithOpNonsingular::mat_mswo_ptr_t
MatrixSymWithOpNonsingular::clone_mswo() const
{
	return clone_mswons();
}

// Overridden from MatrixSymNonsingular

MatrixSymWithOpNonsingular::mat_msns_mut_ptr_t
MatrixSymWithOpNonsingular::clone_msns()
{
	return clone_mswons();
}

MatrixSymWithOpNonsingular::mat_msns_ptr_t
MatrixSymWithOpNonsingular::clone_msns() const
{
	return clone_mswons();
}

// Overridden from MatrixWithOpNonsingular

MatrixSymWithOpNonsingular::mat_mwons_mut_ptr_t
MatrixSymWithOpNonsingular::clone_mwons()
{
	return clone_mswons();
}

MatrixSymWithOpNonsingular::mat_mwons_ptr_t
MatrixSymWithOpNonsingular::clone_mwons() const
{
	return clone_mswons();
}

}	// end namespace AbstractLinAlgPack
