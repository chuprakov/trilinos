// /////////////////////////////////////////////////////////////////////////////////////////////
// MatrixIdentConcat.hpp
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

#ifndef MATRIX_IDENT_CONCAT_H
#define MATRIX_IDENT_CONCAT_H

#include "ConstrainedOptimizationPackTypes.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOp.hpp"

namespace ConstrainedOptimizationPack {

///
/** Matrix class for a matrix vertically concatonated with an identity matrix {abstract}.
 *
 * Represents an interface for a matrix that represents:
 \verbatim

 M = [ alpha*op(D) ]
     [      I      ]
 where:
     D_rng = [1,rows(op(D))]
     I_rng = [rows(op(D))+1,rows(op(D))+cols(op(D))]

 or

 M = [      I      ]
     [ alpha*op(D) ]
 where:
     D_rng = [cols(op(D))+1,rows(op(D))+cols(op(D))]
     I_rng = [1,cols(op(D))]
 \endverbatim
 * and \c I is a <tt>op(D).cols() x op(D).cols()</tt> indentity matrix and
 * the full matrix \c M is of order <tt>(op(D).rows() + op(D).cols()) x op(D).cols()</tt>.
 */
class MatrixIdentConcat
	: virtual public AbstractLinAlgPack::MatrixWithOp
{
public:

	/** @name Access to representation.
	 */
	//@{
	///
	virtual Range1D D_rng() const = 0;
	///
	virtual Range1D I_rng() const = 0;
	///
	virtual value_type alpha() const = 0;
	///
	virtual const MatrixWithOp& D() const = 0;
	///
	virtual BLAS_Cpp::Transp D_trans() const = 0;
	//@}

	/** @name Overridden from MatrixBase */
	//@{
	///
	size_type rows() const;
	///
	size_type cols() const;
	///
	size_type nz() const;
	//@}

	/** @name Overridden from MatrixWithOp */
	//@{
	///
	std::ostream& output(std::ostream& out) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const VectorWithOp& vs_rhs2, value_type beta
		) const;
	///
	void Vp_StMtV(
		VectorWithOpMutable* vs_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const SpVectorSlice& sv_rhs2, value_type beta
		) const;
	//@}

}; // end class MatrixIdentConcat

} // end namespace ConstrainedOptimizationPack

#endif // MATRIX_IDENT_CONCAT_H
