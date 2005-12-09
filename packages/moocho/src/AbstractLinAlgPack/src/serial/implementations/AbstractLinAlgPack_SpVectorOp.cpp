// //////////////////////////////////////////////////////////////////////
// SpVectorOp.cpp
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

#include "AbstractLinAlgPack_SpVectorOp.hpp"

namespace {
// Setup some template classes to check at complile time that
// the layout of SpVectorSlice::element_type is proper.
template<int N, class T>
class assert_compile_time {
public:
assert_compile_time()
{
	// This should not compile if instantiated with a type T that
	// is not an integer.  However, if the compiler checks this
	// function without instantiating it, it can not cause an error
	// because it does not know the type of T to see if the
	// conversion is legal or not.
	T d;
	static_cast<int*>(d);
}
};
// Template specialization for no error
template<>
class assert_compile_time<0,double> {
public:
assert_compile_time()
{}
};
// Validate that there is an integer stride between values
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
	 % (int)sizeof(DenseLinAlgPack::value_type))
	, double
	>
    validate_value_stride;
// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
	 % (int)sizeof(DenseLinAlgPack::index_type))
	, double
	>
    validate_index_stride;
} // end namespace

void AbstractLinAlgPack::add_elements( SpVector* sv_lhs, value_type alpha, const DVectorSlice& vs_rhs
									 , size_type offset, bool add_zeros )
{
	typedef SpVector::element_type ele_t;
	const bool assume_sorted = !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() );
	DVectorSlice::const_iterator
		itr = vs_rhs.begin();
	if(add_zeros) {
		for( size_type i = 1; i <= vs_rhs.dim(); ++i )
			sv_lhs->add_element( ele_t( i + offset, alpha * (*itr++) ) );
	}
	else {
		for( size_type i = 1; i <= vs_rhs.dim(); ++i, ++itr )
			if( *itr != 0.0 )
				sv_lhs->add_element( ele_t( i + offset, alpha * (*itr) ) );
	}
	sv_lhs->assume_sorted(assume_sorted);
}

void AbstractLinAlgPack::add_elements( SpVector* sv_lhs, value_type alpha, const SpVectorSlice& sv_rhs
									 , size_type offset, bool add_zeros )
{
	typedef SpVector::element_type ele_t;
	const bool assume_sorted = ( !sv_lhs->nz() || ( sv_lhs->nz() && sv_lhs->is_sorted() ) )
		&& ( !sv_rhs.nz() || ( sv_rhs.nz() || sv_rhs.is_sorted() ) );
	if(add_zeros) {
		for( SpVectorSlice::const_iterator itr = sv_rhs.begin(); itr != sv_rhs.end(); ++itr )
			sv_lhs->add_element( ele_t( itr->index() + sv_rhs.offset() + offset, alpha * (itr->value()) ) );
	}
	else {
		for( SpVectorSlice::const_iterator itr = sv_rhs.begin(); itr != sv_rhs.end(); ++itr )
			if(itr->value() != 0.0 )
				sv_lhs->add_element( ele_t( itr->index() + sv_rhs.offset() + offset, alpha * (itr->value()) ) );
	}
	sv_lhs->assume_sorted(assume_sorted);
}
