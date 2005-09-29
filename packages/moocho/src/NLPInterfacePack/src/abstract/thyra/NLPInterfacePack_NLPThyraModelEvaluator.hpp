// ///////////////////////////////////////////////////////////////////
// NLPThyraModelEvaluator.hpp
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

#ifndef NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP
#define NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP

#include <vector>

#include "NLPInterfacePack/src/abstract/interfaces/NLPFirstOrder.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorSpace.hpp"
//#include "TSFCoreNonlinTypes.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Teuchos_TestForException.hpp"

namespace NLPInterfacePack {

///
/** Implement the %NLPFirstOrder interface using a
 * <tt>TSFCore::Nonlin::NonlinearProblemFirstOrder</tt> object.
 *
 * This subclass allows great flexibility in how a
 * <tt>TSFCore::Nonlin::NonlinearProblemFirstOrder</tt> object
 * <tt>np</tt> can be used to implement an nonlinear program.
 *
 * The nonlinear program is mapped as follows:

 \verbatim

    min     f(x)
    s.t.    c(x) = 0
            xL <= x <= xu

    where:
            [ np.y                                  ]
        x = [ np.u(u_indep_ind[0])                  ]
            [ ...                                   ]
            [ np.u(u_indep_ind[num_u_indep_sets-1]) ]

        f(x) = sum( q[k], k = 0...num_obj-1 )
        
                                        / g(obj_ind[k])(y,{u(l)})    : if obj_pow[k] == OBJ_LINEAR
            where: q[k] =  obj_wgt[k] * |
                                        \ g(obj_ind[k])(y,{u(l)})^2  : if obj_pow[k] == OBJ_SQUARED

        c(x) = np.c(y,{u(l)})
 \endverbatim

 * where <tt>num_u_indep_sets</tt>, <tt>u_indep_ind[]</tt>,
 * <tt>num_obj</tt>, <tt>obj_ind[]</tt>, <tt>obj_wgt[]</tt> and
 * <tt>obj_pow[]</tt> are input, along with the object <tt>np</tt> to
 * the constructor <tt>NLPThyraModelEvaluator()</tt> or the initialization
 * function <tt>initialize()</tt>.
 *
 * The transformation shown above allows great flexibility is how the
 * auxiliary variables <tt>{u(l)}</tt> auxiliary response functions
 * <tt>np.g(y,{u(l})</tt> are used.  The client can select one or more of
 * the members in the set of auxiliary variables <tt>{u(l)}</tt> to be
 * used as design variables as specified by <tt>num_u_indep_sets</tt>
 * and <tt>u_indep_ind[]</tt>.  Also, the objective function can be
 * built out of a composition of one or more of the auxiliary
 * functions <tt>np.g(y,{u(l)})</tt> as specified by <tt>num_obj</tt>,
 * <tt>obj_ind[]</tt>, <tt>obj_wgt[]</tt> and <tt>obj_pow[]</tt>.
 * This allows one to ignore certain response functions and to build
 * multi-term objective functions with various weights (i.e. as
 * specified <tt>obj_wgt[]</tt>) and perhaps least squares (i.e. as
 * specified by <tt>obj_pow[]</tt>).
 *
 * In addition, the client can also override the bounds on <tt>y</tt>
 * and <tt>{u(l)}</tt> defined in the object <tt>np</tt>.
 *
 * The current implementation of this class does not allow the use of
 * any of the auxiliary functions <tt>np.g(y,{u(l})</tt> as
 * undecomposed equality constraints or extra general inequality
 * constraints.  This type of functionality can be added when it
 * is needed (just ask for it).
 *
 * ToDo: Finish documentation!
 */
class NLPThyraModelEvaluator : virtual public NLPFirstOrder {
public:
	
	///
	enum EObjPow { OBJ_LINEAR, OBJ_SQUARED };

	/// Initialize to uninitialized
	NLPThyraModelEvaluator();

	///
	/** Calls <tt>initialize()</tt>.
	 */
	NLPThyraModelEvaluator(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
		,const int                                                      num_u_indep_sets
		,const int                                                      u_indep_ind[]
		,const int                                                      num_obj
		,const int                                                      obj_ind[]
		,const value_type                                               obj_wgt[]
		,const EObjPow                                                  obj_pow[]
		,const Thyra::VectorBase<value_type>                            *yL       = NULL
		,const Thyra::VectorBase<value_type>                            *yU       = NULL
		,const Thyra::VectorBase<value_type>                            *y0       = NULL
		,const Thyra::VectorBase<value_type>*                           uL[]      = NULL
		,const Thyra::VectorBase<value_type>*                           uU[]      = NULL
		,const Thyra::VectorBase<value_type>*                           u0[]      = NULL
		);

	///
	/** Initialize given a <tt>Thyra::ModelEvaluator</tt> and
	 * a description of how to interpret it.
	 *
	 * @param  np    [in] NonlinearProblem that defines all of the functions.
	 * @param  num_u_indep_sets
	 *               [in] Number of set of auxiliary variables u(l) to include in the set of
	 *               independent variables xI (see above).
	 * @param  num_obj
	 *               [in] The number of auxiliary response functions to include as terms in the
	 *               objective function (see above).
	 * @param  obj_ind
	 *               [in] Array (of length num_obj) that gives the indexes of the auxiliary
	 *               response functions in <tt>g(y,{u(l)})</tt> that will be used as terms
	 *               in the objective function (see above).
	 * @param  obj_wgt
	 *               [in] Array (of length num_obj) that gives the weights for the auxiliary
	 *               response functions in <tt>g(y,{u(l)})</tt> that will be used as terms
	 *               in the objective function (see above).
	 * @param  obj_pow
	 *               [in] Array (of length num_obj) that gives the power to raise for the auxiliary
	 *               response functions in <tt>g(y,{u(l)})</tt> that will be used as terms
	 *               in the objective function (see above).
	 * @param  yL    [in] Pointer to upper bounds for the state variables <tt>y</tt>.  If NULL
	 *               then the default supplied in <tt>np->yL()</tt> will be used.
	 * @param  yU    [in] Pointer to upper bounds for the state variables <tt>y</tt>.  If NULL
	 *               then the default supplied in <tt>np->yL()</tt> will be used.
	 * @param  y0    [in] Pointer to initial guess for the state variables <tt>y</tt>.  If NULL
	 *               the the default supplied in <tt>np->y0()</tt> will be used.
	 * @param  uL    [in] Array (length num_u_indep_sets) of pointers to lower bounds for the
	 *               auxiliary variables <tt>u(l)</tt>.  If NULL then the default supplied in
	 *               <tt>np->uL(u_indep_ind[k])</tt> will be used.
	 * @param  uU    [in] Array (length num_u_indep_sets) of pointers to upper bounds for the
	 *               auxiliary variables <tt>u(l)</tt>.  If NULL then the default supplied in
	 *               <tt>np->uU(u_indep_ind[k])</tt> will be used.
	 * @param  u0    [in] Array (length num_u_indep_sets) of pointers to initial guesses for the
	 *               auxiliary variables <tt>u(l)</tt>.  If NULL then the default supplied in
	 *               <tt>np->u0(u_indep_ind[k])</tt> will be used.
	 *
	 * ToDo: Finish documentation!
	 *
	 * Todo: Add arguments for auxiliary inequalites and equalities
	 */
	void initialize(
		const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
		,const int                                                      num_u_indep_sets
		,const int                                                      u_indep_ind[]
		,const int                                                      num_obj
		,const int                                                      obj_ind[]
		,const value_type                                               obj_wgt[]
		,const EObjPow                                                  obj_pow[]
		,const Thyra::VectorBase<value_type>                            *yL       = NULL
		,const Thyra::VectorBase<value_type>                            *yU       = NULL
		,const Thyra::VectorBase<value_type>                            *y0       = NULL
		,const Thyra::VectorBase<value_type>*                           uL[]      = NULL
		,const Thyra::VectorBase<value_type>*                           uU[]      = NULL
		,const Thyra::VectorBase<value_type>*                           u0[]      = NULL
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
	void imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const;

	//@}

private:

	// /////////////////////////////////////////
	// Private data members

	bool                                initialized_;  // flag for if initialized has been called.
	value_type                          obj_scale_;    // default = 1.0;
	bool                                has_bounds_;   // True if has bounds
	bool                                force_xinit_in_bounds_; // default = true.
	index_type                          num_bounded_x_;
	Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >
	                                    np_;
	int                                 num_u_indep_sets_;
	std::vector<int>                    u_indep_ind_;
	int                                 num_obj_;
	std::vector<int>                    obj_ind_;
	std::vector<value_type>             obj_wgt_;
	std::vector<EObjPow>                obj_pow_;
	std::vector<Range1D>                var_indep_u_;
	bool                                f_has_sqr_term_;
	VectorSpace::space_ptr_t            space_x_;      // Space for the variables
	VectorSpace::space_ptr_t            space_c_;      // Space for the constraints
	NLPFirstOrder::mat_fcty_ptr_t       factory_Gc_;   // Factory for Gc
	NLPFirstOrder::basis_sys_ptr_t      basis_sys_;    // The basis system
	VectorSpace::vec_mut_ptr_t          xinit_;        // Initial guess.
	VectorSpace::vec_mut_ptr_t          xl_;           // lower bounds.
	VectorSpace::vec_mut_ptr_t          xu_;           // upper bounds.

	mutable std::vector<const Thyra::VectorBase<value_type>*>  np_u_in_;

	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >                     np_g_;
	mutable bool                                                             np_g_updated_;
	Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >                     np_c_;
	mutable bool                                                             np_c_updated_;
	Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> >                np_DgDy_;
	std::vector<Teuchos::RefCountPtr<Thyra::MultiVectorBase<value_type> > >  np_DgDu_;
	mutable bool                                                             np_Dg_updated_;
	Teuchos::RefCountPtr<Thyra::LinearOpWithSolveBase<value_type> >          np_DcDy_;
	std::vector<Teuchos::RefCountPtr<Thyra::LinearOpBase<value_type> > >     np_DcDu_;
	mutable bool                                                             np_Dc_updated_;

	mutable bool                                                             f_calc_new_last_;

	// /////////////////////////////////////////
	// Private member functions

	///
	void assert_is_initialized() const;
	///
	void copy_from_y( const Thyra::VectorBase<value_type>& y, VectorMutable* x_D );
	///
	void copy_from_u( const Thyra::VectorBase<value_type> &u, const Range1D& var_indep_u, VectorMutable* x_I );
	///
	void set_x(
		const Vector& x, bool newx
		,const Thyra::VectorBase<value_type>**  y
		,const Thyra::VectorBase<value_type>*** u
		) const;
	///
	void calc_g( const Thyra::VectorBase<value_type> &y, const Thyra::VectorBase<value_type>* u[], bool newx  ) const;
	///
	void calc_Dg( const Thyra::VectorBase<value_type> &y, const Thyra::VectorBase<value_type>* u[], bool newx  ) const;

};	// end class NLPThyraModelEvaluator

// ///////////////////////////////////////////////
// Inline member functions

inline
void NLPThyraModelEvaluator::assert_is_initialized() const
{
	TEST_FOR_EXCEPTION(
		!is_initialized(), NLP::UnInitialized
		,"NLPThyraModelEvaluator::assert_is_initialized() : Error, "
		"NLPThyraModelEvaluator::initialize() has not been called yet."
		);
}

}	// end namespace NLPInterfacePack

#endif	// NLPIP_NLP_THYRA_MODEL_EVALUATOR_HPP
