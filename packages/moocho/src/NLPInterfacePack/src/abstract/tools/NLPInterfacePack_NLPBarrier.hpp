// //////////////////////////////////////////
// BarrierNLP.h
//
// Copyright (C) 2001
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

#ifndef BARRIER_NLP_H
#define BARRIER_NLP_H

#include "NLPInterfacePack/include/NLPObjGradient.h"

namespace NLPInterfacePack {

///
/** Simple wrapper that provides an objective fn with the barrier term
 *   included.
 *
 */
class BarrierNLP
	: public NLPObjGradient
	{
	private:
		MemMngPack::ref_count_ptr<NLPObjGradient> nlp_;
		
	public:
		
		/** @name Public Methods */
		//@{

		/// Set the barrier parameter
		void mu(const value_type mu);

		/// Get the barrier term
		// must be called after calc_f
		value_type barrier_term() const;

		/// Get the true objective term value
		// must be called after calc_f
		value_type objective_term() const;

		/// Get the value of the gradient of the barrier term
		// must be called after calc_Gf
		const MemMngPack::ref_count_ptr<VectorWithOp> grad_barrier_term() const;

		/// Get the value of the gradient of the true objective term
		// must be called after calc_Gf
		const MemMngPack::ref_count_ptr<VectorWithOp> grad_objective_term() const;
		
		//@}
			
		/** @name Constructors / initializers */
		//@{

		///
		/** Constructor.
		 *
		 */
		BarrierNLP();

		///
		void InitializeFromNLP(
		  MemMngPack::ref_count_ptr<NLP> original_nlp
		  );

		// @}

		/** @name Overridden public members from NLPObjGradient */
		//@{

		///
		void initialize(bool test_setup)
		{ nlp_->initialize(test_setup); }
	
		///
		virtual void set_Gf(VectorWithOpMutable* Gf)
		{ nlp_->set_Gf(Gf); }

		///
		virtual VectorWithOpMutable* get_Gf()
		{ return nlp_->get_Gf(); }

		///
		virtual VectorWithOpMutable& Gf()
		{ return nlp_->Gf(); }

		///
		virtual const VectorWithOp& Gf() const
		{ return nlp_->Gf(); }

		/// Overloaded to include barrier term
		virtual void calc_Gf(const VectorWithOp& x, bool newx = true) const;

		///
		virtual size_type num_Gf_evals() const
		{ return nlp_->num_Gf_evals(); }

		//@}

		/** @name Overridden public members from NLP */
		//@{

		///
		//static value_type infinite_bound()

		///
		virtual void force_xinit_in_bounds(bool force_xinit_in_bounds)
		{ nlp_->force_xinit_in_bounds(force_xinit_in_bounds); }

		///
		virtual bool force_xinit_in_bounds() const
		{ return nlp_->force_xinit_in_bounds(); }

		///
		virtual bool is_initialized() const
		{ return nlp_->is_initialized(); }

		///
		virtual size_type n() const
		{ return nlp_->n(); }

		///
		virtual size_type m() const
		{ return nlp_->m(); }

		///
		virtual size_type mI() const
		{ return nlp_->mI(); }
	
		///
		virtual vec_space_ptr_t space_x() const
		{ return nlp_->space_x(); }

		///
		virtual vec_space_ptr_t space_c() const
		{ return nlp_->space_c(); }

		///
		virtual vec_space_ptr_t space_h() const
		{ return nlp_->space_h(); }

		///
		virtual size_type num_bounded_x() const
		{ return nlp_->num_bounded_x(); }

		///
		virtual const VectorWithOp& xl() const
		{ return nlp_->xl(); }

		///
		virtual const VectorWithOp& xu() const 
		{ return nlp_->xu(); }

		///
		virtual value_type max_var_bounds_viol() const
		{ return nlp_->max_var_bounds_viol(); }

		///
		virtual const VectorWithOp& hl() const
		{ return nlp_->hl(); }

		///
		virtual const VectorWithOp& hu() const
		{ return nlp_->hu(); }

		///
		virtual const VectorWithOp& xinit() const
		{ return nlp_->xinit(); }

		///
		virtual void get_init_lagrange_mult(
		  VectorWithOpMutable*   lambda
		  ,VectorWithOpMutable*  lambdaI
		  ,VectorWithOpMutable*  nu
		  ) const
		{ nlp_->get_init_lagrange_mult(lambda, lambdaI, nu); }

		///
		virtual void set_f(value_type* f)
		{ nlp_->set_f(f); }

		///
		virtual value_type* get_f()
		{ return nlp_->get_f(); }

		///
		virtual value_type& f()
		{ return nlp_->f(); }

		///
		virtual const value_type& f() const
		{ return nlp_->f(); }

		///
		virtual void set_c(VectorWithOpMutable* c)
		{ nlp_->set_c(c); }

		///
		virtual VectorWithOpMutable* get_c()
		{ return nlp_->get_c(); }

		///
		virtual VectorWithOpMutable& c()
		{ return nlp_->c(); }

		///
		virtual const VectorWithOp& c() const
		{ return nlp_->c(); }

		///
		virtual void set_h(VectorWithOpMutable* h)
		{ nlp_->set_h(h); }

		///
		virtual VectorWithOpMutable* get_h()
		{ return nlp_->get_h(); }

		///
		virtual VectorWithOpMutable& h()
		{ return nlp_->h(); }

		///
		virtual const VectorWithOp& h() const
		{ return nlp_->h(); }

		///
		virtual void set_multi_calc(bool multi_calc) const
		{ nlp_->set_multi_calc(multi_calc); }

		///
		virtual bool multi_calc() const
		{ return nlp_->multi_calc(); }

		///
		virtual void scale_f( value_type scale_f )
		{ nlp_->scale_f(); }

		///
		virtual value_type scale_f() const
		{ return nlp_->scale_f(); }

		/// Overloaded to include barrier term
		virtual void calc_f(const VectorWithOp& x, bool newx =  true) const;

		///
		virtual void calc_c(const VectorWithOp& x, bool newx = true) const
		{ nlp_->calc_c(x, newx); }

		///
		virtual void calc_h(const VectorWithOp& x, bool newx = true) const
		{ nlp_->calc_h(x, newx); }

		///
		virtual void report_final_solution(
		  const VectorWithOp&    x
		  ,const VectorWithOp*   lambda
		  ,const VectorWithOp*   lambdaI
		  ,const VectorWithOp*   nu
		  ,bool                  is_optimal
		  ) const
		{ nlp_->report_final_solution(
		  x, lambda, lambdaI, nu, is_optimal
		  );
		}

		///
		virtual size_type num_f_evals() const
		{ return nlp_->num_f_evals(); }

		///
		virtual size_type num_c_evals() const
		{ return nlp_->num_c_evals(); }

		///
		virtual size_type num_h_evals() const
		{ return nlp_->num_h_evals(); }

		///
		const ZeroOrderInfo zero_order_info() const
		{ return nlp_->zero_order_info(); }

		//@}

	protected:
		/* protected members Overridden from NLP */
		//@{
		///
		virtual void imp_calc_f(
		  const VectorWithOp& x, 
		  bool newx, 
		  const ZeroOrderInfo& zero_order_info
		  ) const;

		///
		virtual void imp_calc_c(
		  const VectorWithOp& x, 
		  bool newx, 
		  const ZeroOrderInfo& zero_order_info
		  ) const;

		///
		virtual void imp_calc_h(
		  const VectorWithOp& x, 
		  bool newx, 
		  const ZeroOrderInfo& zero_order_info
		  ) const;

		//@}


		/* protected members Overridden from NLPObjGradient */
		//@{

		///
		virtual void imp_calc_Gf(
		  const VectorWithOp& x,
		  bool newx, 
		  const ObjGradInfo& obj_grad_info
		  ) const;

		//@}
	private:

		///
		value_type mu_;
		
		///
		mutable value_type barrier_term_;

		///
		mutable value_type objective_term_;

		/// 
		mutable MemMngPack::ref_count_ptr<VectorWithOpMutable> grad_barrier_term_;
		mutable MemMngPack::ref_count_ptr<VectorWithOpMutable> grad_barrier_term_temp_;
		
		/// 
		mutable MemMngPack::ref_count_ptr<VectorWithOpMutable> grad_objective_term_;

		///
		value_type CalculateBarrierTerm(const VectorWithOp& x) 
const;

	};	// end class BarrierNLP

}	// end namespace NLPInterfacePack

#endif	// BARRIER_NLP_H
