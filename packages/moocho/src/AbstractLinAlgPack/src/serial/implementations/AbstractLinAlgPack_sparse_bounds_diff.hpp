// /////////////////////////////////////////////////////////////////////
// sparse_bounds_diff.h

#ifndef SPARSE_BOUNDS_DIFF_H
#define SPARSE_BOUNDS_DIFF_H

#include "SparseLinAlgPackTypes.h" 

namespace SparseLinAlgPack {

///
/** Take the difference between a spare lower bound vector and a dense vector.
  * 
  * r = alpha * ( sv - v )
  *
  * If sign > 0 then alpha = 1.0.
  * If sign < 0 then alpha = -1.0. 
  * 
  * If uplo == upper, then the nonstored elements in
  * sv are +inf, and if uplo == lower then the nonstored
  * elements are lower.
  */
void imp_sparse_bnd_diff(
	  int						sign
	, const SpVectorSlice		&sv
	, BLAS_Cpp::Uplo			uplo
	, const VectorSlice			&v
	, VectorSlice				*r
	);

}	// end namespace SparseLinAlgPack

#endif  // SPARSE_BOUNDS_DIFF_H
