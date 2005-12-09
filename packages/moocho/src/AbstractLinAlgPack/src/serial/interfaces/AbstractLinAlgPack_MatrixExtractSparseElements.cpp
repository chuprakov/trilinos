// //////////////////////////////////////////////////////////
// MatrixExtractSparseElements.cpp
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

#include <assert.h>

#include "AbstractLinAlgPack_MatrixExtractSparseElements.hpp"
#include "Thyra_Range1D.hpp"

namespace AbstractLinAlgPack {

// Overridden from MatrixConvertToSparseFortranCompatible

index_type
MatrixExtractSparseElements::num_nonzeros(
	EExtractRegion        extract_region
	,EElementUniqueness   element_uniqueness
	) const
{
	index_type dl,du;
	get_dl_du(extract_region,&dl,&du);
	return count_nonzeros(
		element_uniqueness,NULL,NULL,Range1D(1,rows()),Range1D(1,cols()),dl,du);
}

void MatrixExtractSparseElements::coor_extract_nonzeros(
	EExtractRegion                extract_region
	,EElementUniqueness           element_uniqueness
	,const index_type             len_Aval
	,value_type                   Aval[]
	,const index_type             len_Aij
	,index_type                   Arow[]
	,index_type                   Acol[]
	,const index_type             row_offset
	,const index_type             col_offset
	 ) const
{
	int dl,du;
	get_dl_du(extract_region,&dl,&du);
	coor_extract_nonzeros(
		element_uniqueness
		,NULL,NULL,Range1D(1,rows()),Range1D(1,cols()),dl,du
		,1.0
		,len_Aval,Aval,len_Aij,Arow,Acol,row_offset,col_offset);
}

// private

void MatrixExtractSparseElements::get_dl_du(
	EExtractRegion extract_region, index_type* dl, index_type* du
	) const
{
	const size_type
		rows = this->rows(),
		cols = this->cols();
	switch(extract_region) {
		case EXTRACT_FULL_MATRIX:
			*dl = -(rows-1);
			*du = +(cols-1);
			break;
		case EXTRACT_UPPER_TRIANGULAR:
			*dl = 0;
			*du = +(cols-1);
			break;
		case EXTRACT_LOWER_TRIANGULAR:
			*dl = -(rows-1);
			*du = 0;
			break;
		default:
			assert(0);
			break;
	}
}


}	// end namespace AbstractLinAlgPack 
