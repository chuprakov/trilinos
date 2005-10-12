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
#include "AbstractLinAlgPack/src/abstract/thyra/ThyraAccessors.hpp"
#include "AbstractLinAlgPack/src/abstract/thyra/VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack/src/abstract/thyra/VectorMutableThyra.hpp"
#include "AbstractLinAlgPack/src/abstract/thyra/MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/BasisSystemComposite.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack/src/abstract/tools/VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack/src/serial/implementations/MatrixSymPosDefCholFactor.hpp"
#include "Thyra_ExplicitVectorView.hpp"
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
	:initialized_(false),obj_scale_(1.0)
{}

NLPThyraModelEvaluator::NLPThyraModelEvaluator(
	const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
	,const int                                                      num_p_indep_sets
	,const int                                                      p_indep_ind[]
	,const int                                                      num_obj
	,const int                                                      obj_ind[]
	,const value_type                                               obj_wgt[]
	,const EObjPow                                                  obj_pow[]
	,const Thyra::VectorBase<value_type>                            *np_xL
	,const Thyra::VectorBase<value_type>                            *np_xU
	,const Thyra::VectorBase<value_type>                            *np_x0
	,const Thyra::VectorBase<value_type>*                           np_pL[]
	,const Thyra::VectorBase<value_type>*                           np_pU[]
	,const Thyra::VectorBase<value_type>*                           np_p0[]
	)
	:initialized_(false),obj_scale_(1.0)
{
	initialize(np,num_p_indep_sets,p_indep_ind,num_obj,obj_ind,obj_wgt,obj_pow,np_xL,np_xU,np_x0,np_pL,np_pU,np_p0);
}

void NLPThyraModelEvaluator::initialize(
	const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
	,const int                                                      num_p_indep_sets
	,const int                                                      p_indep_ind[]
	,const int                                                      num_obj
	,const int                                                      obj_ind[]
	,const value_type                                               obj_wgt[]
	,const EObjPow                                                  obj_pow[]
	,const Thyra::VectorBase<value_type>                            *np_xL
	,const Thyra::VectorBase<value_type>                            *np_xU
	,const Thyra::VectorBase<value_type>                            *np_x0
	,const Thyra::VectorBase<value_type>*                           np_pL[]
	,const Thyra::VectorBase<value_type>*                           np_pU[]
	,const Thyra::VectorBase<value_type>*                           np_p0[]
	)
{

  TEST_FOR_EXCEPT(num_obj > 1); // ToDo: Make more general!

	using Teuchos::dyn_cast;
	using AbstractLinAlgPack::VectorSpaceThyra;
	using AbstractLinAlgPack::VectorMutableThyra;
	using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef ::Thyra::ModelEvaluatorBase MEB;
	
	initialized_     = false;
	np_g_updated_ = np_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;

  //np_inArgs_ = np->createInArgs();
  np_outArgs_ = np->createOutArgs();
  ::Thyra::ModelEvaluatorBase::DerivativeProperties np_W_properties = np_outArgs_.get_W_properties();
	
#ifdef _DEBUG
	const char msg_err[] = "NLPThyraModelEvaluator::initialize(...): Errror!";
	TEST_FOR_EXCEPTION( np.get() == NULL, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( np_W_properties.supportsAdjoint==false, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( np_W_properties.rank==MEB::DERIV_RANK_DEFICIENT, std::invalid_argument, msg_err );
#endif
	//if(!np->isInitialized()) np->initialize();
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !(num_p_indep_sets <= np->Np()), std::invalid_argument, msg_err );
#endif
	TEST_FOR_EXCEPT( np->Ng() != 1 ); // ToDo: Handle the more general case later!
	const int numResponseFunctions = np->get_g_space(1)->dim(); // ToDo: Handle the more general case later!
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( !(num_obj <= numResponseFunctions), std::invalid_argument, msg_err );
#endif
	//
	//np->initialize(true);
	np_ = np;
	num_p_indep_sets_ = num_p_indep_sets;
	p_indep_ind_.resize(num_p_indep_sets); if(num_p_indep_sets) std::copy( p_indep_ind, p_indep_ind+num_p_indep_sets, p_indep_ind_.begin() );
	num_obj_ = num_obj;
	obj_ind_.resize(num_obj); if(num_obj) std::copy( obj_ind, obj_ind+num_obj, obj_ind_.begin() );
	obj_wgt_.resize(num_obj); if(num_obj) std::copy( obj_wgt, obj_wgt+num_obj, obj_wgt_.begin() );
	obj_pow_.resize(num_obj); if(num_obj) std::copy( obj_pow, obj_pow+num_obj, obj_pow_.begin() );
	//
	f_has_sqr_term_ = false; {for(int p=0; p<num_obj; ++p) if(obj_pow_[p]==OBJ_SQUARED) f_has_sqr_term_=true;}
	
	//
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > np_space_x(np_->get_x_space());
	bool no_y = (np_space_x.get() == NULL);

	//
	std::vector<VectorSpace::space_ptr_t> space_u(num_p_indep_sets);
	if(1){for( int k=0; k<num_p_indep_sets; ++k ) {
		space_u[k] = Teuchos::rcp(new VectorSpaceThyra(np_->get_p_space(p_indep_ind[k])));
	}}
	VectorSpace::space_ptr_t space_xI;
	if(num_p_indep_sets) {
		if(num_p_indep_sets==1) space_xI = space_u[0];
		else                    space_xI = Teuchos::rcp(new VectorSpaceBlocked(&space_u[0],num_p_indep_sets));
	}
	VectorSpace::space_ptr_t space_xD;
	//
	if(!no_y) {
		space_xD = Teuchos::rcp(new VectorSpaceThyra(np_space_x));
		if (num_p_indep_sets)  {
			VectorSpace::space_ptr_t spaces_xD_xI[] = { space_xD, space_xI };
			space_x_ = Teuchos::rcp(new VectorSpaceBlocked(spaces_xD_xI,2));
		}
		else {
			space_x_ = space_xD;
		}
	}
	else {
		space_x_ = space_xI;
	} 

	assert(space_x_.get());

	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > np_space_f(np_->get_f_space());
	bool no_c = (np_space_f.get() == NULL);
	//
	if(!no_c)
		space_c_ = Teuchos::rcp(new VectorSpaceThyra(np_space_f));
	else
		space_c_ = Teuchos::null;

	//
	var_indep_u_.resize(num_p_indep_sets);
	if(1){
		index_type p = 0;
		for( int k=0; k<num_p_indep_sets; ++k ) {
			const index_type space_u_k_dim = space_u[k]->dim();
			var_indep_u_[k] = Range1D(p+1,p+space_u_k_dim);
			p += space_u_k_dim;
		}
	}
	//
	xinit_ = space_x_->create_member();  *xinit_ = 0.0;
	xl_    = space_x_->create_member();  *xl_    = -NLP::infinite_bound();
	xu_    = space_x_->create_member();  *xu_    = +NLP::infinite_bound();

	if(!no_c) {

		factory_Gc_ = BasisSystemComposite::factory_Gc();

		basis_sys_ = Teuchos::rcp(
			new BasisSystemComposite(
				space_x_
				,space_c_
				,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixOpNonsing,MatrixOpNonsingThyra>())          // factory_C
				,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOp,MatrixSymPosDefCholFactor>())         // factory_transDtD
				,Teuchos::rcp(new Teuchos::AbstractFactoryStd<MatrixSymOpNonsing,MatrixSymPosDefCholFactor>())  // factory_S
				)
			);

		if(!no_y) {
			VectorSpace::vec_mut_ptr_t xinit_D = xinit_->sub_view(basis_sys_->var_dep());
			if(np_x0)  copy_from_y( np_x0, xinit_D.get() );
			else       copy_from_y( np_->get_x_init().get(), &*xinit_D );
			VectorSpace::vec_mut_ptr_t xl_D = xl_->sub_view(basis_sys_->var_dep());
			if(np_xL)  copy_from_y( np_xL, xl_D.get() );
			else       copy_from_y( np_->get_x_lower_bounds().get(), &*xl_D );
			VectorSpace::vec_mut_ptr_t xu_D = xu_->sub_view(basis_sys_->var_dep());
			if(np_xU)  copy_from_y( np_xU, xu_D.get() );
			else       copy_from_y( np_->get_x_upper_bounds().get(), &*xu_D );
		}

		if(num_p_indep_sets) {
			VectorSpace::vec_mut_ptr_t xinit_I = xinit_->sub_view(basis_sys_->var_indep());
			if(np_p0) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_p0[k], var_indep_u_[k], &*xinit_I ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_init(p_indep_ind[k]).get(), var_indep_u_[k], &*xinit_I ); }
			VectorSpace::vec_mut_ptr_t xl_I = xl_->sub_view(basis_sys_->var_indep());
			if(np_pL) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_pL[k], var_indep_u_[k], xl_I.get() ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_lower_bounds(p_indep_ind[k]).get(), var_indep_u_[k], &*xl_I ); }
			VectorSpace::vec_mut_ptr_t xu_I = xu_->sub_view(basis_sys_->var_indep());
			if(np_pU) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_pU[k], var_indep_u_[k], xu_I.get() ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_upper_bounds(p_indep_ind[k]).get(), var_indep_u_[k], &*xu_I ); }
		}
	}
  else {

		factory_Gc_ = Teuchos::null;

		basis_sys_ = Teuchos::null;

		if(num_p_indep_sets) {
			VectorSpace::vec_mut_ptr_t xinit_I = xinit_;
			if(np_p0) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_p0[k], var_indep_u_[k], &*xinit_I ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_init(p_indep_ind[k]).get(), var_indep_u_[k], &*xinit_I ); }
			VectorSpace::vec_mut_ptr_t xl_I = xl_;
			if(np_pL) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_pL[k], var_indep_u_[k], xl_I.get() ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_lower_bounds(p_indep_ind[k]).get(), var_indep_u_[k], &*xl_I ); }
			VectorSpace::vec_mut_ptr_t xu_I = xu_;
			if(np_pU) { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_pU[k], var_indep_u_[k], xu_I.get() ); }
			else      { for(int k=0; k<num_p_indep_sets; ++k) copy_from_u( np_->get_p_upper_bounds(p_indep_ind[k]).get(), var_indep_u_[k], &*xu_I ); }
		}

	}
  
	if(num_obj) {
		np_g_ = ::Thyra::createMember(np_->get_g_space(1)); // ToDo: Make more general?
		if(!no_y)
			np_DgDx_ = ::Thyra::createMembers(np_space_x,numResponseFunctions);
		else
			np_DgDx_ = Teuchos::null;
		np_DgDp_.resize(num_p_indep_sets_);
		if(1){for(int k=0; k<num_p_indep_sets; ++k) np_DgDp_[k] = ::Thyra::createMembers(np_->get_p_space(p_indep_ind[k]),numResponseFunctions);}
	}

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
	NLP::set_f(f);
}

void NLPThyraModelEvaluator::set_c(VectorMutable* c)
{
	NLP::set_c(c);
}

void NLPThyraModelEvaluator::unset_quantities()
{
  NLP::unset_quantities();
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
  // ToDo: Do something with this stuff like save it as local data or write it
  // to a file!  Or, add a function to ModelEvaluator that accepts a final
  // point.
}

// Overridden public members from NLPObjGrad

void NLPThyraModelEvaluator::set_Gf(VectorMutable* Gf)
{
  NLPObjGrad::set_Gf(Gf);
}

// Overridden public members from NLPFirstOrder

void NLPThyraModelEvaluator::set_Gc(MatrixOp* Gc)
{
  NLPFirstOrder::set_Gc(Gc);
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
  evalModel(x,newx,&zero_order_info,NULL,NULL);
}

void NLPThyraModelEvaluator::imp_calc_c(
	const Vector& x, bool newx
	,const ZeroOrderInfo& zero_order_info
  ) const
{
  evalModel(x,newx,&zero_order_info,NULL,NULL);
}

// Overridden protected members from NLPObjGrad

void NLPThyraModelEvaluator::imp_calc_Gf(
	const Vector& x, bool newx
	,const ObjGradInfo& obj_grad_info
  ) const
{
  evalModel(x,newx,NULL,&obj_grad_info,NULL);
}

// Overridden protected members from NLPFirstOrder

void NLPThyraModelEvaluator::imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const
{
  evalModel(x,newx,NULL,NULL,&first_order_info);
}

// private

void NLPThyraModelEvaluator::assert_is_initialized() const
{
	TEST_FOR_EXCEPTION(
		!is_initialized(), NLP::UnInitialized
		,"NLPThyraModelEvaluator::assert_is_initialized() : Error, "
		"NLPThyraModelEvaluator::initialize() has not been called yet."
		);
}

void NLPThyraModelEvaluator::copy_from_y( const Thyra::VectorBase<value_type>* y, VectorMutable* x_D ) const
{
  if(!y) return;
	*x_D = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(y),false));
}

void NLPThyraModelEvaluator::copy_from_u( const Thyra::VectorBase<value_type> *u, const Range1D& var_indep_u, VectorMutable* x_I ) const
{
  if(!u) return;
	*x_I->sub_view(var_indep_u) = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(u),false));
}

void NLPThyraModelEvaluator::evalModel( 
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    *zero_order_info
  ,const ObjGradInfo      *obj_grad_info
  ,const FirstOrderInfo   *first_order_info
  ) const
{
  using Teuchos::RefCountPtr;
  using Teuchos::dyn_cast;
  using AbstractLinAlgPack::VectorMutableThyra;
  //
	if(newx) np_g_updated_ = np_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;
  //
  // Set the input arguments
  //
  Thyra::ModelEvaluatorBase::InArgs<value_type>
    np_inArgs = np_->createInArgs();
	if( basis_sys_.get() != NULL) {
		const Range1D
			var_dep   = basis_sys_->var_dep(),
			var_indep = basis_sys_->var_indep();
		RefCountPtr<const Vector> xD = x.sub_view(var_dep), xI;
		if(num_p_indep_sets_) xI = x.sub_view(var_indep);
    np_inArgs.set_x(dyn_cast<const VectorMutableThyra>(*xD).thyra_vec());
		if(num_p_indep_sets_) {
			for( int k=0; k<num_p_indep_sets_; ++k ) {
        np_inArgs.set_p(p_indep_ind_[k],dyn_cast<const VectorMutableThyra>(*xI->sub_view(var_indep_u_[k])).thyra_vec());
			}
    }
	}
	else { // no dependent vars
		assert(num_p_indep_sets_);
		for( int k=0; k<num_p_indep_sets_; ++k ) {
      np_inArgs.set_p(p_indep_ind_[k],dyn_cast<const VectorMutableThyra>(*x.sub_view(var_indep_u_[k])).thyra_vec());
		}
	}
  //
  // Set the output arguments
  //
  MatrixOp            *Gc = NULL;
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  //
  if(zero_order_info) {
    f = zero_order_info->f;
    c = zero_order_info->c;
  }
  else if(obj_grad_info) {
    Gf = obj_grad_info->Gf;
    f  = obj_grad_info->f;
    c  = obj_grad_info->c;
  }
  else if(first_order_info) {
    Gc = first_order_info->Gc;
    Gf = first_order_info->Gf;
    f  = first_order_info->f;
    c  = first_order_info->c;
  }
  else {
    TEST_FOR_EXCEPT(true); // Should never be called!
  }
  //
  Thyra::ModelEvaluatorBase::OutArgs<value_type> np_outArgs = np_->createOutArgs();
  if( f && num_obj_ && !f_updated_ ) {
    np_outArgs.set_g(1,np_g_); // ToDo: Make more general!
  }
  if( c && !c_updated_ ) {
    Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > thyra_c;
    get_thyra_vector(*space_c_,c,&thyra_c);
    np_outArgs.set_f(thyra_c);
  }
  if( Gf && !Gf_updated_ ) {
    TEST_FOR_EXCEPT(true);
  }
  if( Gc && !Gc_updated_ ) {
    TEST_FOR_EXCEPT(true);
  }
  //
  // Evaluate the model
  //
  np_->evalModel(np_inArgs,np_outArgs);
  //
  // Update outputs
  //
  if( f && !f_updated_ ) {
    if(num_obj_) {
      *f = ::Thyra::get_ele(*np_g_,1);
      np_outArgs.set_g(1,Teuchos::null);
    }
    else {
      *f = 0.0;
    }
    f_updated_ = true;
  }
  if( c && !c_updated_ ) {
    Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >
      thyra_c = np_outArgs.get_f();
    commit_thyra_vector(*space_c_,c,&thyra_c);
    np_outArgs.set_f(Teuchos::null);
    c_updated_ = true;
  }
  if( Gf && !Gf_updated_ ) {
    TEST_FOR_EXCEPT(true);
  }
  if( Gc && !Gc_updated_ ) {
    TEST_FOR_EXCEPT(true);
  }
}

}	// end namespace NLPInterfacePack
