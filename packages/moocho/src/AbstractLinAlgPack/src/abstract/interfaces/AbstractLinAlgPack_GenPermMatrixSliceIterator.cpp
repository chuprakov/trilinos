// /////////////////////////////////////////////////////////
// GenPermMatrixSliceIterator.cpp
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

#include <string>

#include "AbstractLinAlgPack_GenPermMatrixSliceIterator.hpp"
#include "Teuchos_TestForException.hpp"

void AbstractLinAlgPack::GenPermMatrixSliceIteratorPack::GPMS_row_col_iterator_assert_not_null(
	const void* p)
{
	TEST_FOR_EXCEPTION(
		!p, std::logic_error
		,"GenPermMatrixSliceIteratorPack::row_col_iterator<T>, Error "
		"row_i can not be NULL" );
}
