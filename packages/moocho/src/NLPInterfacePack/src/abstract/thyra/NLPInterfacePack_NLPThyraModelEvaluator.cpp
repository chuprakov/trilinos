// ///////////////////////////////////////////////////////////////////
// NLPThyraModelEvaluator.cpp
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

#include <assert.h>

#include <algorithm>

#include "NLPThyraModelEvaluator.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/LinAlgOpPack.hpp"
#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
//#include "AbstractLinAlgPack/src/abstract/tsfcore/VectorSpaceTSFCore.hpp"
//#include "AbstractLinAlgPack/src/abstract/tsfcore/VectorMutableTSFCore.hpp"
//#include "AbstractLinAlgPack/src/abstract/tsfcore/MatrixOpNonsingTSFCore.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/BasisSystemComposite.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/serial/implementations/MatrixSymPosDefCholFactor.hpp"
#include "TSFCoreNonlinLinearOpWithSolve.hpp"
#include "TSFCoreNonlinNonlinearProblemFirstOrder.hpp"
#include "TSFCoreExplicitVectorView.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

// Debugging only
//#include <iostream>
//#include <typeinfo>
//#include "AbstractLinAlgPack/src/abstract/interfaces/VectorOut.hpp"
//#include "TSFCoreTestingTools.hpp"

namespace NLPInterfacePack {

NLPThyraModelEvaluator::NLPThyraModelEvaluator()
	:initialized_(false),obj_scale_(1.0),f_calc_new_last_(false)
{}

NLPThyraModelEvaluator::NLPThyraModelEvaluator(
	const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
	,const int                                                      num_u_indep_sets
	,const int                                                      u_indep_ind[]
	,const int                                                      num_obj
	,const int                                                      obj_ind[]
	,const value_type                                               obj_wgt[]
	,const EObjPow                                                  obj_pow[]
	,const Thyra::VectorBase<value_type>                            *yL
	,const Thyra::VectorBase<value_type>                            *yU
	,const Thyra::VectorBase<value_type>                            *y0
	,const Thyra::VectorBase<value_type>*                           uL[]
	,const Thyra::VectorBase<value_type>*                           uU[]
	,const Thyra::VectorBase<value_type>*                           u0[]
	)
	:initialized_(false),obj_scale_(1.0),f_calc_new_last_(false)
{
	initialize(np,num_u_indep_sets,u_indep_ind,num_obj,obj_ind,obj_wgt,obj_pow,yL,yU,y0,uL,uU,u0);
}

/**  If we have no dependent variables y, we set basis_sys_ = null;
 */

void NLPThyraModelEvaluator::initialize(
	const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
	,const int                                                      num_u_indep_sets
	,const int                                                      u_indep_ind[]
	,const int                                                      num_obj
	,const int                                                      obj_ind[]
	,const value_type                                               obj_wgt[]
	,const EObjPow                                                  obj_pow[]
	,const Thyra::VectorBase<value_type>                            *yL
	,const Thyra::VectorBase<value_type>                            *yU
	,const Thyra::VectorBase<value_type>                            *y0
	,const Thyra::VectorBase<value_type>*                           uL[]
	,const Thyra::VectorBase<value_type>*                           uU[]
	,const Thyra::VectorBase<value_type>*                           u0[]
	)
{
  TEST_FOR_EXCEPT(0);
}
	
// Overridden public members from NLP

void NLPThyraModelEvaluator::initialize(bool test_setup)
{
	if(initialized_) {
		NLPFirstOrder::initialize(test_setup);
		return;
	}
	//assert(0); // Todo: push the variables in bounds!
	num_bounded_x_ = AbstractLinAlgPack::num_bounded(*xl_,*xu_,NLP::infinite_bound());
	NLPFirstOrder::initialize(test_setup);
	initialized_ = true;
}

bool NLPThyraModelEvaluator::is_initialized() const
{
	return initialized_;
}

NLP::vec_space_ptr_t
NLPThyraModelEvaluator::space_x() const
{
	return space_x_;
}

NLP::vec_space_ptr_t
NLPThyraModelEvaluator::space_c() const
{
	return space_c_;
}

size_type NLPThyraModelEvaluator::num_bounded_x() const
{
	return num_bounded_x_;
}

void NLPThyraModelEvaluator::force_xinit_in_bounds(bool force_xinit_in_bounds)
{
	force_xinit_in_bounds_ = force_xinit_in_bounds;
}

bool NLPThyraModelEvaluator::force_xinit_in_bounds() const
{
	return force_xinit_in_bounds_;
}

const Vector& NLPThyraModelEvaluator::xinit() const
{
	return *xinit_;
}

const Vector& NLPThyraModelEvaluator::xl() const
{
	return *xl_;
}

const Vector& NLPThyraModelEvaluator::xu() const
{
	return *xu_;
}

value_type NLPThyraModelEvaluator::max_var_bounds_viol() const
{
	return 1e-5; // I have no idea?
}

void NLPThyraModelEvaluator::set_f(value_type* f)
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::set_c(VectorMutable* c)
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::unset_quantities()
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::scale_f( value_type scale_f )
{
	obj_scale_ = scale_f;
}

value_type NLPThyraModelEvaluator::scale_f() const
{
	return obj_scale_;
}

void NLPThyraModelEvaluator::report_final_solution(
	const Vector&    x
	,const Vector*   lambda
	,const Vector*   nu
	,bool            optimal
	)
{
  TEST_FOR_EXCEPT(0);
}

// Overridden public members from NLPObjGrad

void NLPThyraModelEvaluator::set_Gf(VectorMutable* Gf)
{
  TEST_FOR_EXCEPT(0);
}

// Overridden public members from NLPFirstOrder

void NLPThyraModelEvaluator::set_Gc(MatrixOp* Gc)
{
  TEST_FOR_EXCEPT(0);
}

const NLPFirstOrder::mat_fcty_ptr_t
NLPThyraModelEvaluator::factory_Gc() const
{
	return factory_Gc_;
}

const NLPFirstOrder::basis_sys_ptr_t
NLPThyraModelEvaluator::basis_sys() const
{
	return basis_sys_;
}

// Overridden protected members from NLP

void NLPThyraModelEvaluator::imp_calc_f(
	const Vector& x, bool newx
	,const ZeroOrderInfo& zero_order_info
  ) const
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::imp_calc_c(
	const Vector& x, bool newx
	,const ZeroOrderInfo& zero_order_info
  ) const
{
  TEST_FOR_EXCEPT(0);
}

// Overridden protected members from NLPObjGrad

void NLPThyraModelEvaluator::imp_calc_Gf(
	const Vector& x, bool newx
	,const ObjGradInfo& obj_grad_info
  ) const
{
  TEST_FOR_EXCEPT(0);
}

// Overridden protected members from NLPFirstOrder

void NLPThyraModelEvaluator::imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const
{
  TEST_FOR_EXCEPT(0);
}

// private

void NLPThyraModelEvaluator::copy_from_y( const Thyra::VectorBase<value_type>& y, VectorMutable* x_D )
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::copy_from_u( const Thyra::VectorBase<value_type> &u, const Range1D& var_indep_u, VectorMutable* x_I )
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::set_x(
	const Vector& x, bool newx
	,const Thyra::VectorBase<value_type>** y
	,const Thyra::VectorBase<value_type>*** u
	) const
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::calc_g(
	const Thyra::VectorBase<value_type>   &y
	,const Thyra::VectorBase<value_type>* u[]
	,bool newx 
	) const
{
  TEST_FOR_EXCEPT(0);
}

void NLPThyraModelEvaluator::calc_Dg(
	const Thyra::VectorBase<value_type>   &y
	,const Thyra::VectorBase<value_type>* u[]
	,bool newx 
	) const
{
  TEST_FOR_EXCEPT(0);
}

}	// end namespace NLPInterfacePack
