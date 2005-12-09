// //////////////////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixBasisNonbasisStd.hpp
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

#ifndef MATRIX_BASIS_NON_BASIS_STD_H
#define MATRIX_BASIS_NON_BASIS_STD_H

#include "AbstractLinAlgPack_MatrixBasisNonbasis.hpp"
#include "Miref_count_ptr.h"

namespace AbstractLinAlgPack {

///
/** This is the interface to a basis and an nonbasis matrix.
  *
  * The form of this matrix is:
  *
  * M = [ C' ; N' ]
  */
class MatrixBasisNonbasisStd
	: public MatrixBasisNonbasis
{
public:
	
    ///
	typedef Teuchos::RefCountPtr<const MatrixOp>
		C_ptr_t;
    ///
	typedef Teuchos::RefCountPtr<const MatrixFactorized>
		C_nonsingular_ptr_t;
	///
	typedef Teuchos::RefCountPtr<const MatrixOp>
		N_ptr_t;

	///
	/** Initialize to null pointers for C, C_nonsingular and N.
	 *
	 * Postconditions:\begin{itemize}
	 * \item this->rows() == 0
	 * \item this->cols() == 0
	 * \end{itemize}
	 *
	 * The operations this->C() and this->N() should not be called
	 * until this->initialize(...) has been called with non
	 * -null arguments successfully.
	 */
	MatrixBasisNonbasisStd();

	///
	/** Initialize smart reference counted points to C, C_nonsingular and N.
	  *
	  * Preconditions:\begin{itemize}
	  * \item [C.get() != NULL] C->rows() == C->cols() (throw std::length_error)
	  * \item [C_nonsingular.get() != NULL]
	  *                          C_nonsingular->rows() == C->rows()
	  *                          && C_nonsingular->cols() == C->cols() (throw std::length_error)
	  * \item [N.get() != NULL] C->rows() == N->rows() (throw std::length_error)
	  * \item [N.get() != NULL] C.get() != NULL (throw std::logic_error)
	  * \item [C_nonsingular.get() != NULL] C.get() != NULL (throw std::logic_error)
	  * \end{itemize}
	  *
	  * Postconditions:\begin{itemize}
	  * \item this->rows() == C->cols() + N->cols()
	  * \item this->cols() == C->rows()
	  * \item [C_nonsingular.get() == NULL] this->C_nonsingular() throws std::logic_error
	  * \end{itemize}
      *
	  * Note that it is legal for C.get()!=NULL and N.get()==NULL and therefore
	  * *this = [ C' ] only.
	  *
	  * The behavior of the rest of this object could be spelled out but
	  * the properties should be obvious.  The matrix *this acts like the
	  * concatenated matrix [ C' ; N' ]
	  *
	  * To uninitialize, set C.get()==NULL, N.get()==NULL.
	  *
	  * Of course since ref_count_prt<...> objects are being used for
	  * C and N, calling initialize(...) again or deleting *this object
	  * will cause the current C and N matrices to be deleted if their
	  * reference counts are zero.
	  *
      */
	void initialize( const C_ptr_t& C, const C_nonsingular_ptr_t& C_nonsingular
	    , const N_ptr_t& N );

	///
	C_ptr_t C_ptr() const {
		return C_;
	}
	///
	C_nonsingular_ptr_t C_nonsingular_ptr() const {
		return C_nonsingular_;
	}
	///
	N_ptr_t N_ptr() const {
		return N_;
	}

	// /////////////////////////////////////////////////////
	// Overridden from Matrix

	///
	size_type rows() const;

	///
	size_type cols() const;

	// //////////////////////////////////////////
	// Overridden from MatrixBasisNonbasis

	///
	const MatrixOp& C() const;
	///
	const MatrixFactorized& C_nonsingular() const;
	///
	const MatrixOp& N() const;

private:
	C_ptr_t               C_;
	C_nonsingular_ptr_t   C_nonsingular_;
	N_ptr_t               N_;

	// Validate that the matrix is initialized (throws std::logic_error)
	void assert_initialized() const;

};	// end class MatrixBasisNonbasisStd

}	// end namespace AbstractLinAlgPack 

#endif	// MATRIX_BASIS_NON_BASIS_STD_H
