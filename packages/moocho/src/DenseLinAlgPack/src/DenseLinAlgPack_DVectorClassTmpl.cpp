// //////////////////////////////////////////////////////////////////////////////////
// VectorClassTmpl.cpp

#include <iomanip>

#include "LinAlgPack/include/VectorClassTmpl.h"

#ifdef LINALGPACK_CHECK_SLICE_SETUP
LinAlgPack::size_type LinAlgPack::vector_validate_sized(size_type size)
{
	if(!size)
		throw std::invalid_argument("vector_validate_sized(...) : Error, A vector region can not be created from an unsized vector.\n");
	return size;
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void LinAlgPack::vector_validate_range(size_type ubound, size_type max_ubound)
{
	if(ubound > max_ubound)
		throw std::out_of_range("vector_validate_range(...) : The upper bound is out of range.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RANGE
void LinAlgPack::vector_validate_subscript(size_type size, size_type i)
{
	if(i < 1 || i > size) 
		throw std::out_of_range("vector_validate_subscript(size,i) : Error, Subscript i out of bounds.\n");
}
#endif

#ifdef LINALGPACK_CHECK_RHS_SIZES
void LinAlgPack::assert_vs_sizes(size_type size1, size_type size2)
{
	if(size1 != size2)
		throw std::length_error("assert_vs_sizes(...) : Error, Sizes of vector regions must be the same.\n");
}
#endif
