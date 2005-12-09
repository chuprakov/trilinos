// ////////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixScaling_Strategy.hpp
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

#ifndef MATRIX_SCALING_STRATEGY_H
#define MATRIX_SCALING_STRATEGY_H

namespace AbstractLinAlgPack {

///
/** Abstract interface for sparse matrix scaling strategies
 *
 * ToDo: Finish documentation!
 */
class MatrixScaling_Strategy {
public:

	///
	virtual ~MatrixScaling_Strategy() {}

	/** @name Pure virtual methods to be overridden by subclasses */
	//@{

	/// Scale the matrix and save the scalings for later use for rhs and lhs.
	virtual void scale_matrix(
		index_type m, index_type n, index_type nz
		,const index_type row_i[], const index_type col_j[]
		,bool new_matrix, value_type A[]
		) = 0;
	
	/// Scale the rhs vector
	virtual void scale_rhs( BLAS_Cpp::Transp trans, value_type b[] ) const = 0;
	
	/// Scale the lhs vector
	virtual void scale_lhs( BLAS_Cpp::Transp trans, value_type x[] ) const = 0;
	
	//@}

}; // end class MatrixScaling_Strategy

} // end namespace AbstractLinAlgPack

#endif // MATRIX_SCALING_STRATEGY_H
