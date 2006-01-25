// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP

#include <vector>

#include "NLPInterfacePack_NLPFirstOrder.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
//#include "TSFCoreNonlinTypes.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_TestForException.hpp"

namespace AbstractLinAlgPack { class VectorSpaceThyra; }

namespace NLPInterfacePack {

/** \brief Implement the %NLPFirstOrder interface using a
 * <tt>Thyra::ModelEvaluator</tt> object.
 *
 * The nonlinear program is mapped as follows:

 \verbatim

    min     f(x)
    s.t.    c(x) = 0
            xL <= x <= xu

    where:

        x = [ model.x        ]
            [ model.p(p_idx) ]

        f(x) = model.g(g_idx)

        c(x) = model.f()

 \endverbatim

 * where <tt>p_idx > 0</tt> and <tt>g_idx > 0</tt> are some indexes that
 * specificy the indepenent variables and the objective function.
 *
 * In addition, the client can also override the bounds on <tt>model.x</tt> and
 * <tt>model.p(p_idx)</tt> defined in the object <tt>model</tt>.
 *
 * The current implementation of this class does not allow the use of any of
 * the auxiliary functions <tt>model.g()</tt> as undecomposed equality
 * constraints or extra general inequality constraints.  This type of
 * functionality can be added when it is needed (just ask for it).
 *
 * ToDo: Finish documentation!
 */
class NLPThyraModelEvaluator : virtual public NLPFirstOrder {
public:
	
	///
	enum EDerivSupport { DERIV_SUPPORT_ADJOINT_ONLY, DERIV_SUPPORT_ADJOINT_OR_DIRECT };

	/// Initialize to uninitialized
	NLPThyraModelEvaluator();

	///
	/** Calls <tt>initialize()</tt>.
	 */
	NLPThyraModelEvaluator(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
		,const Thyra::VectorBase<value_type>                            *model_xL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_xU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_x0      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_p0      = NULL
		);

	///
	/** Initialize given a <tt>Thyra::ModelEvaluator</tt> and
	 * a description of how to interpret it.
	 *
	 * @param  model    [in] NonlinearProblem that defines all of the functions and variables.
   * @param  p_idx [in] Index of the subset of parameter vectors to use as the independent
   *               variables.  If <tt>p_idx < 0</tt>, then no extra parameters are added.
   * @param  g_idx [in] Index of the subset of auxiliary response functions to use as
   *               the objective function.  Note, only the first element <tt>model.g(g_idx)(1)</tt>
   *               will be used as the objective function value.
	 * @param  model_xL [in] Pointer to upper bounds for the state variables <tt>model.x</tt>.  If NULL
	 *               then the default supplied in <tt>model->get_x_lower_bounds()</tt> will be used.
	 * @param  model_xU [in] Pointer to upper bounds for the state variables <tt>x</tt>.  If NULL
	 *               then the default supplied in <tt>model->get_x_upper_bounds()</tt> will be used.
	 * @param  model_x0 [in] Pointer to initial guess for the state variables <tt>x</tt>.  If NULL
	 *               the the default supplied in <tt>model->get_x_init()</tt> will be used.
	 *
	 * ToDo: Finish documentation!
	 *
	 * Todo: Add arguments for auxiliary inequalites and equalities
	 */
	void initialize(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &model
    ,const int                                                      p_idx
    ,const int                                                      g_idx
		,const Thyra::VectorBase<value_type>                            *model_xL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_xU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_x0      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pL      = NULL
		,const Thyra::VectorBase<value_type>                            *model_pU      = NULL
		,const Thyra::VectorBase<value_type>                            *model_p0      = NULL
		);

	/** @name Overridden public members from NLP */
	//@{

	///
	void initialize(bool test_setup);
	///
	bool is_initialized() const;
	///
	vec_space_ptr_t space_x() const;
	///
	vec_space_ptr_t space_c() const;
	///
    size_type num_bounded_x() const;
	///
	void force_xinit_in_bounds(bool force_xinit_in_bounds);
	///
	bool force_xinit_in_bounds() const;
	///
	const Vector& xinit() const;
	///
	const Vector& xl() const;
	///
	const Vector& xu() const;
	///
	value_type max_var_bounds_viol() const;
	///
	void set_f(value_type* f);
	///
	void set_c(VectorMutable* c);
	///
	void unset_quantities();
	///
	void scale_f( value_type scale_f );
	///
	value_type scale_f() const;
	///
	void report_final_solution(
		const Vector&    x
		,const Vector*   lambda
		,const Vector*   nu
		,bool            optimal
		);

	//@}

	/** @name Overridden public members from NLPObjGrad */
	//@{

	///
	void set_Gf(VectorMutable* Gf);

	//@}

	/** @name Overridden public members from NLPFirstOrder */
	//@{

	/// Overridden to check the concrete type of Gc
	void set_Gc(MatrixOp* Gc);
	///
	const NLPFirstOrder::mat_fcty_ptr_t factory_Gc() const;
	/// Returns an ExampleBasisSystem
	const basis_sys_ptr_t basis_sys() const;

	//@}

protected:

	/** @name Overridden protected members from NLP */
	//@{

	///
	void imp_calc_f(
		const Vector& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;
	///
	void imp_calc_c(
		const Vector& x, bool newx
		,const ZeroOrderInfo& zero_order_info) const;

	//@}

	/** @name Overridden protected members from NLPObjGrad */
	//@{

	///
	void imp_calc_Gf(
		const Vector& x, bool newx
		,const ObjGradInfo& obj_grad_info) const;

	//@}

	/** @name Overridden protected members from NLPFirstOrder */
	//@{

	///
	void imp_calc_Gc(
    const Vector& x, bool newx
    ,const FirstOrderInfo& first_order_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private types

  typedef Teuchos::RefCountPtr<const AbstractLinAlgPack::VectorSpaceThyra> VectorSpaceThyra_ptr_t;

	// /////////////////////////////////////////
	// Private data members

	bool                                initialized_;  // flag for if initialized has been called.
	value_type                          obj_scale_;    // default = 1.0;
	bool                                has_bounds_;   // True if has bounds
	bool                                force_xinit_in_bounds_; // default = true.
	index_type                          num_bounded_x_;
	Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >
	                                    model_;
  int                                 p_idx_;
  int                                 g_idx_;
  bool                                DfDp_supports_op_;
  bool                                DfDp_supports_mv_;
	VectorSpace::space_ptr_t            space_x_;      // Space for the variables
	VectorSpaceThyra_ptr_t              space_c_;      // Space for the constraints
	NLPFirstOrder::mat_fcty_ptr_t       factory_Gc_;   // Factory for Gc
	NLPFirstOrder::basis_sys_ptr_t      basis_sys_;    // The basis system
	VectorSpace::vec_mut_ptr_t          xinit_;        // Initial guess.
	VectorSpace::vec_mut_ptr_t          xl_;           // lower bounds.
	VectorSpace::vec_mut_ptr_t          xu_;           // upper bounds.

	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >                     model_g_;

  mutable bool model_g_updated_;
  mutable bool model_Dg_updated_;

  mutable bool f_updated_;
  mutable bool c_updated_;
  mutable bool Gf_updated_;
  mutable bool Gc_updated_;

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;
	///
	void copy_from_model_x( const Thyra::VectorBase<value_type>* model_x, VectorMutable* x_D ) const;
	///
	void copy_from_model_p( const Thyra::VectorBase<value_type> *model_p, VectorMutable* x_I ) const;
  ///
  void evalModel( 
		const Vector            &x
    ,bool                   newx
		,const ZeroOrderInfo    *zero_order_info  // != NULL if only zero-order info
    ,const ObjGradInfo      *obj_grad_info    // != NULL if obj-grad and below info
    ,const FirstOrderInfo   *first_order_info // != NULL if first-order and below info
    ) const;

};	// end class NLPThyraModelEvaluator

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP
