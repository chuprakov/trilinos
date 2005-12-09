// /////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_sparse_bounds_diff.hpp
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

#ifndef SPARSE_BOUNDS_DIFF_H
#define SPARSE_BOUNDS_DIFF_H

#include "AbstractLinAlgPack_Types.hpp" 

namespace AbstractLinAlgPack {

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
	, const DVectorSlice			&v
	, DVectorSlice				*r
	);

}	// end namespace AbstractLinAlgPack

#endif  // SPARSE_BOUNDS_DIFF_H
