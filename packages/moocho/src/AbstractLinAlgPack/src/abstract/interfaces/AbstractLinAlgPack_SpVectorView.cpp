// ////////////////////////////////////////////////////////////////////
// SpVectorView.cpp
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

#include "AbstractLinAlgPack/include/SpVectorView.h"

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
	 % (int)sizeof(AbstractLinAlgPack::value_type))
	, double
	>
validate_value_stride;

// Validate that there is an integer stride between indexes
assert_compile_time<
    ((int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
	 % (int)sizeof(AbstractLinAlgPack::index_type))
	, double
	>
validate_index_stride;

// Compute the stride between values
const int values_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
	/ (int)sizeof(AbstractLinAlgPack::value_type);

// Compute the stride between indexes
const int indices_stride = (int)sizeof(AbstractLinAlgPack::SpVectorSlice::element_type)
	/ (int)sizeof(AbstractLinAlgPack::index_type);

} // end namespace

RTOp_SubVector AbstractLinAlgPack::sub_vec_view( const SpVectorSlice& sv )
{
	RTOp_SubVector sub_vec;
	if(!sv.nz()) {
		RTOp_sub_vector_sparse(
			0           // global_offset
			,sv.dim()   // sub_dim
			,0          // nz
			,NULL       // vlaues
			,1          // values_stride
			,NULL       // indices
			,1          // indices_stride
			,0          // local_offset
			,1          // is_sorted
			,&sub_vec);
	}
	else {
		SpVectorSlice::const_iterator
			itr = sv.begin();
		assert(itr != sv.end());
		const value_type   *values  = &itr->value();
		const index_type   *indexes = &itr->index();
		if( sv.nz() == sv.dim() && sv.is_sorted() ) {
			RTOp_sub_vector_dense(
				0                 // global_offset
				,sv.dim()         // sub_dim
				,values           // values
				,values_stride    // values_stride
				,&sub_vec);
		}
		else {
			RTOp_sub_vector_sparse(
				0                 // global_offset
				,sv.dim()         // sub_dim
				,sv.nz()          // nz
				,values           // values
				,values_stride    // values_stride
				,indexes          // indices
				,indices_stride   // indices_stride
				,sv.offset()      // local_offset
				,sv.is_sorted()   // is_sorted
				,&sub_vec);
		}
	}
	
	return sub_vec;
}
