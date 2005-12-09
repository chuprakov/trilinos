// ////////////////////////////////////////////////////////
// ConstrainedOptPack_MatrixHessianSuperBasicInitDiagonal.hpp
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

#ifndef MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H
#define MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H

#include <vector>

#include "ConstrainedOptPack_MatrixHessianSuperBasic.hpp"
#include "AbstractLinAlgPack/src/AbstractLinAlgPack_MatrixSymInitDiag.hpp"

namespace ConstrainedOptPack {

///
/** Matrix class that adds the ability to initialize to a diagonal
 * to a MatrixHessainSuperBasic object.
 *
 * Essentially, the matrix #B_RR# must support the #MatrixSymInitDiag#
 * interface.
 */

class MatrixHessianSuperBasicInitDiagonal
	: public virtual MatrixHessianSuperBasic
	, public virtual MatrixSymInitDiag
	{
public:

	///
	/** Constructs to uninitialized.
	 */
	MatrixHessianSuperBasicInitDiagonal();

	/** @name Overridden from MatrixHessianSuperBasic.
	 **/
	//@{

	///
	/** Initialize the matrix and require B_RR to support #MatrixSymInitDiag#.
	 *
	 * Preconditions:\begin{itemize}
	 * \item #dynamic_cast<MatrixSymInitDiag*>(const_cast<MatrixSymWithOpFactorized*>(B_RR_ptr.get())) != NULL#
	 * \end{itemize}
	 *
	 * This overridden function is a little bit of a hack but it is better than the
	 * alternatives that I could think of.  What is important is that the #MatrixSymInitDiag#
	 * interface for this class will be supported as long as this function executes
	 * succesfully on initialization.  This is far better than waiting until later when
	 * dynamic casting would be done and failed.  We need to catch any mistakes right away.
	 *
	 * In this manner we will have given up some compile-time checking but the run-time check will
	 * be very good.
	 */
	void initialize(
		size_type            n
		,size_type           n_R
		,const size_type     i_x_free[]
		,const size_type     i_x_fixed[]
		,const EBounds       bnd_fixed[]
		,const B_RR_ptr_t&   B_RR_ptr
		,const B_RX_ptr_t&   B_RX_ptr
		,BLAS_Cpp::Transp    B_RX_trans
		,const B_XX_ptr_t&   B_XX_ptr
		);
	
	//@}

	/** @name Overridden from MatrixSymInitDiag.
	 *
	 * These function call the corresponding functions on
	 * #B_RR_ptr()# and then call the equivalent of:
	 \begin{verbatim}
	 MatrixHessianSuperBasic::initialize(
	     this->B_RR_ptr()->rows(),this->B_RR_ptr()->cols()
		 ,NULL,NULL,NULL
		 ,this->B_RR_ptr()
		 ,this->B_RX_ptr(),this->B_RX_trans()
		 ,this->B_XX_ptr()
		);
	\end{verbatim}
	 **/
	//@{

	///
	void init_identity( size_type n, value_type alpha );
	///
	void init_diagonal( const DVectorSlice& diag );

	//@}

private:

	// ///////////////////////////////////
	// Private data members

	MatrixSymInitDiag   *B_RR_init_;  // will be non-null for any valid initalization

	// //////////////////////////
	// Private member functions
	
}; // end class MatrixHessianSuperBasicInitDiagonal

} // end namespace ConstrainedOptPack

#endif // MATRIX_HESSIAN_SUPER_BASIC_INIT_DIAGONAL_H
