// //////////////////////////////////////////////////////////////////////////////////
// COOMatrixTmplConvertToSparseCompressedColumnDecl.h

#ifndef COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
#define COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H

#include "SparseLinAlgPackTypes.h"

namespace SparseLinAlgPack {

/** @name {\bf Conversion to Fortran compatable sparse compressed column 
  * operations for COOMatrixTemplateInterface (Level 2,3 BLAS)}.
  *
  * See the ConvertToSparseCompressedColumn class.
  */
//@{

///
template<class T_COOM>
size_type COOM_num_in_column(
	  const T_COOM&						m
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	);

///
template<class T_COOM>
void COOM_insert_nonzeros(
	  const T_COOM&						m
	, BLAS_Cpp::Transp					trans
	, value_type						alpha
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			);

///
template<class T_COOM>
value_type COOM_insert_scaled_nonzeros(
	  const T_COOM&						m
	, BLAS_Cpp::Transp					trans
	, value_type						scaled_max_ele
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			);

//@}

} // end namespace SparseLinAlgPack

#endif	// COO_MATRIX_TMPL_CONVERT_TO_SPARSE_COMPRESSED_COLUMN_DECL_H
