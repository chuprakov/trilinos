// //////////////////////////////////////////
// BarrierNLP.cpp
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

#include <math.h>
#include <iostream>
#include <limits>

#include "NLPInterfacePack/src/BarrierNLP.hpp"
#include "AbstractLinAlgPack/src/VectorSpace.hpp"
#include "AbstractLinAlgPack/src/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/VectorWithOpOut.hpp"
#include "ThrowException.hpp"

namespace NLPInterfacePack {

BarrierNLP::BarrierNLP()
	:
	barrier_term_(0.0),
	objective_term_(0.0),
	nlp_(MemMngPack::null)
	{
	}


void BarrierNLP::InitializeFromNLP(
  MemMngPack::ref_count_ptr<NLP> original_nlp
  )
	{
	THROW_EXCEPTION(
	  !original_nlp.get(),
	  std::logic_error,
	  "null nlp passed to BarrierNLP decorator"
	  );

	nlp_ = MemMngPack::rcp_dynamic_cast<NLPObjGradient>(original_nlp);

	THROW_EXCEPTION(
	  !nlp_.get(),
	  std::logic_error,
	  "non NLPObjGradient NLP passed to BarrierNLP decorator"
	  );
	}

void BarrierNLP::mu(const value_type mu)
	{
	mu_ = mu;
	}

value_type BarrierNLP::barrier_term() const
	{
	return barrier_term_;
	}

value_type BarrierNLP::objective_term() const
	{
	return objective_term_;
	}

const MemMngPack::ref_count_ptr<VectorWithOp> BarrierNLP::grad_barrier_term() const
	{
	return grad_barrier_term_;
	}

const MemMngPack::ref_count_ptr<VectorWithOp>  BarrierNLP::grad_objective_term() const
	{
	return grad_objective_term_;
	}


void BarrierNLP::calc_f(const VectorWithOp& x, bool newx) const
	{
	nlp_->calc_f(x, newx);
	value_type* f_val = nlp_->get_f();

	objective_term_ = *f_val;
	barrier_term_   = CalculateBarrierTerm(x);

	(*f_val) += barrier_term_;
	}

void BarrierNLP::calc_Gf(const VectorWithOp& x, bool newx) const
	{
	using AbstractLinAlgPack::inv_of_difference;

   	nlp_->calc_Gf(x, newx);
	grad_objective_term_ = nlp_->get_Gf()->clone();

	//std::cout << "grad_objective_term=\n";
	//grad_objective_term_->output(std::cout);

	if (!grad_barrier_term_temp_.get())
		{ grad_barrier_term_temp_ = grad_objective_term_->space().create_member(); }

	if (!grad_barrier_term_.get())
		{ grad_barrier_term_ = grad_objective_term_->space().create_member(); }	

	*grad_barrier_term_temp_ = 0.0;
	*grad_barrier_term_ = 0.0;	

	inv_of_difference(mu_, nlp_->xu(), x, grad_barrier_term_.get());
 	//std::cout << "mu*invXU=\n";
	//grad_barrier_term_->output(std::cout);

	inv_of_difference(mu_, x, nlp_->xl(), grad_barrier_term_temp_.get());
 	//std::cout << "mu*invXL=\n";
	//grad_barrier_term_temp_->output(std::cout);

	grad_barrier_term_->axpy(-1.0, *grad_barrier_term_temp_);

	nlp_->get_Gf()->axpy(1.0, *grad_barrier_term_);

	//std::cout << "grad_objective_term with barrier=\n";
	//nlp_->get_Gf()->output(std::cout);
	}

void BarrierNLP::imp_calc_f(
  const VectorWithOp& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
	{
	assert(false && !"This should never get called.");
	}

void BarrierNLP::imp_calc_c(
  const VectorWithOp& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
	{
	assert(false && !"This should never get called.");
	}

void BarrierNLP::imp_calc_h(
  const VectorWithOp& x, 
  bool newx, 
  const ZeroOrderInfo& zero_order_info
  ) const
	{	
	assert(false && !"This should never get called.");
	}

void BarrierNLP::imp_calc_Gf(
  const VectorWithOp& x,
  bool newx, 
  const ObjGradInfo& obj_grad_info
  ) const
	{
	assert(false && !"This should never get called.");
	}


value_type BarrierNLP::CalculateBarrierTerm(const VectorWithOp& x) const
	{
	using AbstractLinAlgPack::log_bound_barrier;
	barrier_term_ = log_bound_barrier(x, xl(), xu());
//	std::cerr << "BarrierNLP::CalculateBarrierTerm(x) : (1) barrier_term_ = " << barrier_term_ << std::endl;
	barrier_term_ *= -mu_;
//	std::cerr << "BarrierNLP::CalculateBarrierTerm(x) : (2) barrier_term_ = " << barrier_term_ << std::endl;
	return barrier_term_;
	}

}	// end namespace NLPInterfacePack
