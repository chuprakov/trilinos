// /////////////////////////////////////////////////////////////////////
// sparse_bounds_diff.cpp

#include "SparseLinAlgPack/include/sparse_bounds_diff.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/LinAlgOpPack.h"
#include "LinAlgPack/include/LinAlgPackAssertOp.h"

void SparseLinAlgPack::imp_sparse_bnd_diff(
	  int						sign
	, const SpVectorSlice		&sv
	, BLAS_Cpp::Uplo			uplo
	, const VectorSlice			&v
	, VectorSlice				*r
	)
{
	LinAlgPack::Vp_V_assert_sizes(r->size(),sv.size());
	LinAlgPack::VopV_assert_sizes(sv.size(),v.size());

	typedef LinAlgPack::value_type value_type;
	const value_type
		inf = std::numeric_limits<value_type>::max();
	*r = ( uplo == BLAS_Cpp::upper ? inf : -inf );
	const SparseLinAlgPack::SpVectorSlice::difference_type o = sv.offset();
	for( SparseLinAlgPack::SpVectorSlice::const_iterator itr = sv.begin();
			itr != sv.end(); ++itr )
	{
		(*r)(itr->indice() + o) = itr->value();
	}
	LinAlgPack::Vp_StV( r, -1.0, v );
	if( sign < 0 )
		LinAlgPack::Vt_S( r, -1.0 );
}
