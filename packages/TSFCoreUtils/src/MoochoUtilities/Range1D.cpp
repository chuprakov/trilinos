// ////////////////////////////////////////////////////////////////////////////
// Range1D.cpp
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

#include "Range1D.hpp"
#include "Teuchos_TestForException.hpp"

namespace RangePack {

const Range1D Range1D::Invalid(Range1D::INVALID);

//#ifdef _DEBUG

void Range1D::assert_valid_range(size_t lbound, size_t ubound) const {
	TEST_FOR_EXCEPTION(
		lbound < 1, std::range_error
		,"Range1D::assert_valid_range(): Error, lbound ="<<lbound<<" must be greater than 0." );
	TEST_FOR_EXCEPTION(
		lbound > ubound, std::range_error
		,"Range1D::assert_valid_range(): Error, lbound = "<<lbound<<" > ubound = "<<ubound );
}

//#endif // _DEBUG

} // end namespace RangePack
