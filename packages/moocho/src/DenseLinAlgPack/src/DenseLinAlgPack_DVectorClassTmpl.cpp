// //////////////////////////////////////////////////////////////////////////////////
// DVectorClassTmpl.cpp
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

#include <iomanip>

#include "DenseLinAlgPack_DVectorClassTmpl.hpp"
#include "Teuchos_TestForException.hpp"

#ifdef LINALGPACK_CHECK_SLICE_SETUP
DenseLinAlgPack::size_type DenseLinAlgPack::vector_validate_sized(size_type size)
{
	TEST_FOR_EXCEPTION(
		!size, std::invalid_argument
		,"vector_validate_sized(...) : Error, A vector region can not be created from an unsized vector.\n"
		);
	return size;
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_range(size_type ubound, size_type max_ubound)
{
	TEST_FOR_EXCEPTION(
		ubound > max_ubound, std::out_of_range
		,"vector_validate_range(...) : The upper bound is out of range.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void DenseLinAlgPack::vector_validate_subscript(size_type size, size_type i)
{
	TEST_FOR_EXCEPTION(
		i < 1 || i > size, std::out_of_range
		,"vector_validate_subscript(size,i) : Error, Subscript i out of bounds.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RHS_SIZES
void DenseLinAlgPack::assert_vs_sizes(size_type size1, size_type size2)
{
	TEST_FOR_EXCEPTION(
		size1 != size2, std::length_error
		,"assert_vs_sizes(...) : Error, size1 = " << size1 << " != size2 = " << size2 << "\n");
}
#endif
