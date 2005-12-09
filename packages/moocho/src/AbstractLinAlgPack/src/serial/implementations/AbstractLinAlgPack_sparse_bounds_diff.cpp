// /////////////////////////////////////////////////////////////////////
// sparse_bounds_diff.cpp
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

#include "AbstractLinAlgPack_sparse_bounds_diff.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_SpVectorClass.hpp"
#include "AbstractLinAlgPack_LinAlgOpPackHack.hpp"
#include "DenseLinAlgPack_AssertOp.hpp"

void AbstractLinAlgPack::imp_sparse_bnd_diff(
	  int						sign
	, const SpVectorSlice		&sv
	, BLAS_Cpp::Uplo			uplo
	, const DVectorSlice			&v
	, DVectorSlice				*r
	)
{
	DenseLinAlgPack::Vp_V_assert_sizes(r->size(),sv.size());
	DenseLinAlgPack::VopV_assert_sizes(sv.size(),v.size());

	typedef DenseLinAlgPack::value_type value_type;
	const value_type
		inf = std::numeric_limits<value_type>::max();
	*r = ( uplo == BLAS_Cpp::upper ? inf : -inf );
	const AbstractLinAlgPack::SpVectorSlice::difference_type o = sv.offset();
	for( AbstractLinAlgPack::SpVectorSlice::const_iterator itr = sv.begin();
			itr != sv.end(); ++itr )
	{
		(*r)(itr->indice() + o) = itr->value();
	}
	DenseLinAlgPack::Vp_StV( r, -1.0, v );
	if( sign < 0 )
		DenseLinAlgPack::Vt_S( r, -1.0 );
}
