// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_MatrixOpThyra.hpp
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

#ifndef ALAP_MATRIX_OP_Thyra_HPP
#define ALAP_MATRIX_OP_Thyra_HPP

#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "Thyra_LinearOpBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>MatrixOp</tt> adapter subclass for <tt>Thyra::LinearOpBase</tt>.
 */
class MatrixOpThyra : virtual public MatrixOp {
public:

	/** @name Constructors / Initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_linear_op().get() == NULL</tt>
	 * <li><tt>this->rows() == 0</tt>
	 * <li><tt>this->cols() == 0</tt>
	 * </ul>
	 */
	MatrixOpThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MatrixOpThyra(
		const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
		,BLAS_Cpp::Transp                                                    thyra_linear_op_trans = BLAS_Cpp::no_trans
		);
	///
	/** Initalize given a smart pointer to a <tt>Thyra::LinearOpBase</tt> object.
	 *
	 * @param  thyra_linear_op  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_linear_op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li><tt>thyra_linear_op->opSupported(Thyra::NOTRANS) && thyra_linear_op->opSupported(Thyra::TRANS)
	 *     (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_linear_op().get() == thyra_linear_op.get()</tt>
	 * <li><tt>this->thyra_linear_op_trans() == thyra_linear_op_trans</tt>
	 * <li><tt>this->rows() == thyra_linear_op->range()->dim()</tt>
	 * <li><tt>this->cols() == thyra_linear_op->domain()->dim()</tt>
	 * </ul>
	 */
	virtual void initialize(
		const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >   &thyra_linear_op
		,BLAS_Cpp::Transp                                                    thyra_linear_op_trans = BLAS_Cpp::no_trans
		);
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::VectorBase</tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_linear_op().get() == NULL</tt>
	 * </ul>
	 */
	Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> > set_uninitialized();
	///
	/** Return a (converted) smart pointer to the internal smart pointer to the <tt>Thyra::VectorBase</tt> object.
	 *
	 * If <tt>this->thyra_linear_op().count() == 1</tt>, then <tt>this</tt>
	 * has sole ownership of the <tt>*this->thyra_linear_op()</tt> object.
	 */
	const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >& thyra_linear_op() const;
	///
	BLAS_Cpp::Transp thyra_linear_op_trans() const;

	//@}

	/** @name Overridden from MatrixBase */
	//@{

	///
	const VectorSpace& space_cols() const;
	///
	const VectorSpace& space_rows() const;

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	///
	mat_mut_ptr_t clone();
	///
	MatrixOp& operator=(const MatrixOp& mwo_rhs);
	///
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2, value_type beta
		) const;
	/// Works for <tt>MultiVectorMutableThyra</tt> arguments
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta
		) const;

	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	Thyra::LinearOpBase<value_type>                                *thyra_linear_op;
#else
	Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >   thyra_linear_op_;
#endif
	BLAS_Cpp::Transp                                               thyra_linear_op_trans_;
	VectorSpaceThyra                                               space_cols_;
	VectorSpaceThyra                                               space_rows_;

}; // end class MatrixOpThyra

// //////////////////////////////////////////////////
// Inlined functions

inline
const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >&
MatrixOpThyra::thyra_linear_op() const
{
	return thyra_linear_op_;
}

inline
BLAS_Cpp::Transp MatrixOpThyra::thyra_linear_op_trans() const
{
	return thyra_linear_op_trans_;
}

} // end namespace AbstractLinAlgPack

#endif // ALAP_MATRIX_OP_Thyra_HPP
