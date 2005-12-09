// ////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_ConvertToCSC.hpp
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

#ifndef CONVERT_TO_SPARSE_COMPRESSED_COLUMN_H
#define CONVERT_TO_SPARSE_COMPRESSED_COLUMN_H

#include "AbstractLinAlgPack_Types.hpp"
#include "DenseLinAlgPack_IVector.hpp"

/** @name Conversion utilites / interfaces for a Fortran compatable sparse compressed column matrix format.
  *
  * These declarations and functions are ment to support the creation
  * of a Fortran compatable sparse matrix in compressed column
  * format given as set of matrix objects of arbirary format.\\
  *
    \begin{verbatim}
  	        [               a2*op(A2)       ]
  	D = P * [    a1*op(A1)                  ] * Q
  	        [               ai*op(Ai)       ]
    \end{verbatim}
  *
  * Where:\begin{itemize}
  *		\item	Ai are as set of matrix objects (dense vectors or matrices, 
  *				or ConvertToCSC types),
  *		\item	ai are the scaling parameters.
  *		\item	op(Ai) = Ai or Ai' and nonoverlapping regions of the nonpermuted
  *				formed matrix
  *		\item	P and Q are permutation matrices.
  * \end{itemize}
  *
  * These are used to form the compressed column data structure for D
  * composed of the real array D_val(nz), and integer arrays D_col_start(n)
  * and D_row_i(nz) where:
  *
  *	n :			number of columns in D (indexed by k)\\\\
  *	nz :		number of nonzero elements for the entire matrix D\\\\
  *	D_val(k) :	The value of the kth nonzero element sorted by column.\\\\
  *	D_col_start(j):	The location in D_val for the start of the elements
  *				in the jth column of D.  D_col_start(n+1) gives the posstion
  *				in D_val one past the last element in the last column.  This
  *				is so that you can calcuate the number of nonzero elements
  *				in a column j as ( D_col_start(j+1) - D_col_start(j) ).\\\\
  * D_row_i(k):	row indice for the kth nonzero element stored by column as
  *				in D_val.\\
  *
  * Note that since it is allowed to transpose Ai that you can obtain the
  * compressed row representation by applying this method to the D' given
  * above.
  *
  * ToDo: finish documentation.
  */
//@{

namespace AbstractLinAlgPack {

///
/** Abstract interface for inserting a matrix into
  * a Fortran compatable compressed column sparse
  * matrix format.
  */
class ConvertToCSC {
public:

	///
	virtual ~ConvertToCSC()
	{}

	///
	/** Count the number of nonzero elements in each
	  * column of the matrix by incrementing num_in_col
	  *
	  * @return The number of nonzeros added.
	  */
	virtual size_type num_in_column(
		  BLAS_Cpp::Transp					trans
		, size_type							col_offset
		, const IVector::value_type*		col_perm
		, size_type*						num_in_col	) const = 0;
		
	///
	/** Inserts the nozero elements.
	  */
	virtual void insert_nonzeros(
		  BLAS_Cpp::Transp					trans
		, value_type						alpha
		, size_type							row_offset
		, size_type							col_offset
		, const IVector::value_type*		row_perm
		, const IVector::value_type*		col_perm
		, size_type*						next_nz_in_col
		, FortranTypes::f_dbl_prec*			D_val
		, FortranTypes::f_int*				D_row_i			) const = 0;

	///
	/** Inserts the nozero elements scaled so that max|alpha*A| = scaled_max_ele
	  * and returns alpha.
	  */
	virtual value_type insert_scaled_nonzeros(
		  BLAS_Cpp::Transp					trans
		, value_type						scaled_max_ele
		, size_type							row_offset
		, size_type							col_offset
		, const IVector::value_type*		row_perm
		, const IVector::value_type*		col_perm
		, size_type*						next_nz_in_col
		, FortranTypes::f_dbl_prec*			D_val
		, FortranTypes::f_int*				D_row_i			) const = 0;

};	// end class ConvertToCSC

namespace ConvertToSparseCompressedColumnPack {

/** @name Conversion to spare compressed column matrix format for concrete types.
  */
//@{

///
/** Add a nonzero element for a scalar.
  */
inline void scalar_insert_nonzero(
	  value_type						alpha
	, size_type							row_i
	, size_type							col_j
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			)
{
	size_type ele = next_nz_in_col[ col_j - 1 ]++;
	D_val[ ele - 1 ] = alpha;
	if( D_row_i )
		D_row_i[ ele - 1 ] = row_i;
}

///
/** Add the nonzero elements in a dense vector for each column.
  */
void vector_insert_nonzeros(
	  const DVectorSlice&				vs
	, value_type						alpha
	, size_type							row_offset
	, size_type							col_j
	, const IVector::value_type*		row_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			);

///
/** Count the number of nonzeros in a dense matrix for each column.
  */
size_type dense_num_in_column(
	  size_type							rows
	, size_type							cols
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	);

///
/** Add the nonzero elements in a dense matrix for each column.
  */
void dense_insert_nonzeros(
	  const DMatrixSlice&				gms
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
/** Add the nonzero elements in a dense matrix for each column scaled.
  */
value_type dense_insert_scaled_nonzeros(
	  const DMatrixSlice&				gms
	, BLAS_Cpp::Transp					trans
	, value_type						scaled_max_ele
	, size_type							row_offset
	, size_type							col_offset
	, const IVector::value_type*		row_perm
	, const IVector::value_type*		col_perm
	, size_type*						next_nz_in_col
	, FortranTypes::f_dbl_prec*			D_val
	, FortranTypes::f_int*				D_row_i			);

/** @name Add a MatrixOp object.  If the matrix supports the conversion inteface
  * then all is dandy.  If not then a conversion to dense must be preformed.
  */
//@{

///
size_type num_in_column(
	  const MatrixOp&				m
	, BLAS_Cpp::Transp					trans
	, size_type							col_offset
	, const IVector::value_type*		col_perm
	, size_type*						num_in_col	);

///
void insert_nonzeros(
	  const MatrixOp&				m
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
value_type insert_scaled_nonzeros(
	  const MatrixOp&				m
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

//@}

}	// end namespace ConvertToSparseCompressedColumnPack

}	// end namespace AbstractLinAlgPack

//	end Conversion utilities
//@}

// /////////////////////////////////////////////////////////////////
// Inline definitions

#endif	// CONVERT_TO_SPARSE_COMPRESSED_COLUMN_H
