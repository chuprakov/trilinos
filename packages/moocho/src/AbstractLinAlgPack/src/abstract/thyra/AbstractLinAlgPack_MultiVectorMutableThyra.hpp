// ///////////////////////////////////////////////////////////
// AbstractLinAlgPack_MultiVectorMutableThyra.hpp
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

#ifndef ALAP_MULTI_VECTOR_MUTABLE_Thyra_HPP
#define ALAP_MULTI_VECTOR_MUTABLE_Thyra_HPP

#include "AbstractLinAlgPack_MultiVectorMutable.hpp"
#include "AbstractLinAlgPack_MatrixOpThyra.hpp"
#include "Thyra_MultiVectorBase.hpp"

namespace AbstractLinAlgPack {

///
/** <tt>MultiVectorMutable</tt> adapter subclass for <tt>Thyra::MultiVectorBase</tt>.
 */
class MultiVectorMutableThyra
	:virtual public MultiVectorMutable
	,virtual public MatrixOpThyra
{
public:

	/** @name Constructors / Initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_vec().get() == NULL</tt>
	 * <li>See the post conditions for <tt>MatrixOpThyra::MatrixOpThyra()</tt>
	 * </ul>
	 */
	MultiVectorMutableThyra();
	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	MultiVectorMutableThyra( const Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> >& thyra_multi_vec );
	///
	/** Initalize given a smart pointer to a <tt>Thyra::MultiVectorBase</tt> object.
	 *
	 * @param  thyra_multi_vec  [in] Smart pointer to Thyra vector <tt>this</tt> will adapt.
	 *
	 * Preconditioins:<ul>
	 * <li><tt>thyra_multi_vec.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_multi_vec().get() == thyra_multi_vec.get()</tt>
	 * <li>See the post conditions for <tt>MatrixOpThyra::initialize()</tt>
	 * </ul>
	 */
	void initialize( const Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> >& thyra_multi_vec );
	///
	/** Set to uninitialized and return smart pointer to the internal <tt>Thyra::LinearOpBase</tt> object.
	 *
	 * Postconditioins:<ul>
	 * <li><tt>this->thyra_multi_vec().get() == NULL</tt>
	 * </ul>
	 *
	 * Note that his nonvirtual function hides the nonvirtual function
	 * <tt>MatrixOpThyra::set_uninitialized()</tt>.
	 */
	Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> > set_uninitialized();
	///
	/** Return a smart pointer to the internal <tt>Thyra::LinearOpBase</tt> object.
	 */
	Teuchos::RefCountPtr<const Thyra::MultiVectorBase<value_type> > thyra_multi_vec() const;

	//@}

	/** @name Overridden from MatrixOpThyra */
	//@{

	/// Performs a const_cast<> and dynamic_cast<> and passes on to <tt>this->initialize()</tt>.
	void initialize( const Teuchos::RefCountPtr<const Thyra::LinearOpBase<value_type> >& thyra_linear_op );

	//@}

	/** @name Overridden from MatrixOp */
	//@{

	/// Overridden to call <tt>MatrixOpThyra::clone()</tt>
	mat_mut_ptr_t clone();
	/// Overridden to call <tt>MultiVectorMutable::operator=()</tt>
	MatrixOp& operator=(const MatrixOp& mwo_rhs);
	/// Overridden to call <tt>MatrixOpThyra::Vp_StMtV()</tt>
	void Vp_StMtV(
		VectorMutable* v_lhs, value_type alpha, BLAS_Cpp::Transp trans_rhs1
		,const Vector& v_rhs2, value_type beta
		) const;
	/// Overridden to call <tt>MatrixOpThyra::Mp_StMtM()</tt>
	bool Mp_StMtM(
		MatrixOp* mwo_lhs, value_type alpha
		,BLAS_Cpp::Transp trans_rhs1
		,const MatrixOp& mwo_rhs2, BLAS_Cpp::Transp trans_rhs2
		,value_type beta
		) const;

	//@}

	/** @name Overridden from MultiVector */
	//@{

	///  Returns <tt>COL_ACCESS</tt>
	access_by_t access_by() const;
	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op
		,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
		,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
		,RTOpPack::ReductTarget* reduct_objs[]
		,const index_type primary_first_ele,   const index_type primary_sub_dim, const index_type primary_global_offset
		,const index_type secondary_first_ele, const index_type secondary_sub_dim
		) const;
	///
	void apply_op(
		EApplyBy apply_by, const RTOpPack::RTOp& primary_op, const RTOpPack::RTOp& secondary_op
		,const size_t num_multi_vecs,      const MultiVector*   multi_vecs[]
		,const size_t num_targ_multi_vecs, MultiVectorMutable*  targ_multi_vecs[]
		,RTOpPack::ReductTarget *reduct_obj
		,const index_type primary_first_ele, const index_type primary_sub_dim, const index_type primary_global_offset
		,const index_type secondary_first_ele, const index_type secondary_sub_dim
		) const;

	//@}

	/** @name Overridden from MultiVectorMutable */
	//@{

	///
	vec_mut_ptr_t col(index_type j);
	/// <tt>return.get()==NULL</tt>
	vec_mut_ptr_t row(index_type i);
	/// <tt>return.get()==NULL</tt>
	vec_mut_ptr_t diag(int k);
	///
	multi_vec_mut_ptr_t mv_sub_view(const Range1D& row_rng, const Range1D& col_rng);

	//@}

private:

	///
	Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> > cast_thyra_multi_vec();
	
}; // end class MultiVectorMutableThyra

} // end namespace AbstractLinAlgPack

#endif // ALAP_MULTI_VECTOR_MUTABLE_Thyra_HPP
