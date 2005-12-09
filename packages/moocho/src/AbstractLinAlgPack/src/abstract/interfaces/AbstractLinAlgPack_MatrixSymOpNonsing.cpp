// //////////////////////////////////////////////////////////////////////
// MatrixSymOpNonsing.cpp
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

#include "AbstractLinAlgPack_MatrixSymOpNonsing.hpp"

namespace AbstractLinAlgPack {

MatrixSymOpNonsing::mat_mswons_mut_ptr_t
MatrixSymOpNonsing::clone_mswons()
{
	return Teuchos::null;
}

MatrixSymOpNonsing::mat_mswons_ptr_t
MatrixSymOpNonsing::clone_mswons() const
{
	return Teuchos::null;
}

// Overridden from MatrixOp

MatrixSymOpNonsing::mat_mut_ptr_t
MatrixSymOpNonsing::clone()
{
	return clone_mswons();
}

MatrixSymOpNonsing::mat_ptr_t
MatrixSymOpNonsing::clone() const
{
	return clone_mswons();
}

// Overridden from MatrixNonsing

MatrixSymOpNonsing::mat_mns_mut_ptr_t
MatrixSymOpNonsing::clone_mns()
{
	return clone_mswons();
}

MatrixSymOpNonsing::mat_mns_ptr_t
MatrixSymOpNonsing::clone_mns() const
{
	return clone_mswons();
}

// Overridden from MatrixSymOp

MatrixSymOpNonsing::mat_mswo_mut_ptr_t
MatrixSymOpNonsing::clone_mswo()
{
	return clone_mswons();
}

MatrixSymOpNonsing::mat_mswo_ptr_t
MatrixSymOpNonsing::clone_mswo() const
{
	return clone_mswons();
}

// Overridden from MatrixSymNonsing

MatrixSymOpNonsing::mat_msns_mut_ptr_t
MatrixSymOpNonsing::clone_msns()
{
	return clone_mswons();
}

MatrixSymOpNonsing::mat_msns_ptr_t
MatrixSymOpNonsing::clone_msns() const
{
	return clone_mswons();
}

// Overridden from MatrixOpNonsing

MatrixSymOpNonsing::mat_mwons_mut_ptr_t
MatrixSymOpNonsing::clone_mwons()
{
	return clone_mswons();
}

MatrixSymOpNonsing::mat_mwons_ptr_t
MatrixSymOpNonsing::clone_mwons() const
{
	return clone_mswons();
}

}	// end namespace AbstractLinAlgPack
