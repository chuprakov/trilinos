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

//
// Note: I am using a BasisSystem here just to make it easy to generate
// var_dep() and var_indep() and for not much else.
//

#include "NLPInterfacePack_NLPThyraModelEvaluator.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_BasisSystemComposite.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_MatrixSymPosDefCholFactor.hpp"
#include "Thyra_ExplicitVectorView.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {
	
// Overridden public members from NLP

void NLPThyraModelEvaluator::initialize(bool test_setup)
{
	if(initialized_) {
		NLPObjGrad::initialize(test_setup);
		return;
	}
	//assert(0); // Todo: push the variables in bounds!
	num_bounded_x_ = AbstractLinAlgPack::num_bounded(*xl_,*xu_,NLP::infinite_bound());
	NLPObjGrad::initialize(test_setup);
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
  f_updated_ = false;
}

void NLPThyraModelEvaluator::set_c(VectorMutable* c)
{
	NLP::set_c(c);
  c_updated_ = false;
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
  Gf_updated_ = false;
}

// Overridden protected members from NLP

void NLPThyraModelEvaluator::imp_calc_f(
	const Vector& x, bool newx
	,const ZeroOrderInfo& zero_order_info
  ) const
{
  evalModel(x,newx,&zero_order_info,NULL);
}

void NLPThyraModelEvaluator::imp_calc_c(
	const Vector& x, bool newx
	,const ZeroOrderInfo& zero_order_info
  ) const
{
  evalModel(x,newx,&zero_order_info,NULL);
}

// Overridden protected members from NLPObjGrad

void NLPThyraModelEvaluator::imp_calc_Gf(
	const Vector& x, bool newx
	,const ObjGradInfo& obj_grad_info
  ) const
{
  evalModel(x,newx,NULL,&obj_grad_info);
}


// Protected functions to be used by subclasses

NLPThyraModelEvaluator::NLPThyraModelEvaluator()
	:initialized_(false),obj_scale_(1.0)
{}

void NLPThyraModelEvaluator::initializeBase(
  const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &model
  ,const int                                                      p_idx
  ,const int                                                      g_idx
  ,const Thyra::VectorBase<value_type>                            *model_xL
  ,const Thyra::VectorBase<value_type>                            *model_xU
  ,const Thyra::VectorBase<value_type>                            *model_x0
  ,const Thyra::VectorBase<value_type>                            *model_pL
  ,const Thyra::VectorBase<value_type>                            *model_pU
  ,const Thyra::VectorBase<value_type>                            *model_p0
	)
{

	using Teuchos::dyn_cast;
	using AbstractLinAlgPack::VectorSpaceThyra;
	using AbstractLinAlgPack::VectorMutableThyra;
	using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef ::Thyra::ModelEvaluatorBase MEB;
	
	initialized_ = false;
	model_g_updated_ = model_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;

	const char msg_err[] = "NLPThyraModelEvaluator::initialize(...): Errror!";
  Thyra::ModelEvaluatorBase::OutArgs<double> model_outArgs = model->createOutArgs();
	TEST_FOR_EXCEPTION( model.get() == NULL, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( p_idx >= 0 && ( p_idx > model_outArgs.Np()-1 ), std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( g_idx >= 0 && ( g_idx > model_outArgs.Ng()-1 ), std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( !model_outArgs.supports(MEB::OUT_ARG_W), std::invalid_argument, msg_err );
  MEB::DerivativeProperties model_W_properties = model_outArgs.get_W_properties();
	TEST_FOR_EXCEPTION( model_W_properties.supportsAdjoint==false, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( model_W_properties.rank==MEB::DERIV_RANK_DEFICIENT, std::invalid_argument, msg_err );
  if(p_idx >= 0 ) {
    TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).none(), std::invalid_argument, msg_err );
    if(g_idx >= 0) {
      TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DgDp,g_idx,p_idx).none(), std::invalid_argument, msg_err );
    }
  }
  if(g_idx >= 0) {
    TEST_FOR_EXCEPTION( model_outArgs.supports(MEB::OUT_ARG_DgDx,g_idx).none(), std::invalid_argument, msg_err );
  }
	//
	model_ = model;
  p_idx_ = p_idx;
  g_idx_ = g_idx;
  //
  if(p_idx >= 0 ) {
    DfDp_supports_op_ = model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_LINEAR_OP);
    DfDp_supports_mv_ = model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_MV_BY_COL);
  }
  else {
    DfDp_supports_op_ = false;
    DfDp_supports_mv_ = false;
  }
	//
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > model_space_x(model_->get_x_space());
	bool no_model_x = (model_space_x.get() == NULL);
	//
	VectorSpace::space_ptr_t space_xI;
  if(p_idx >= 0)
    space_xI = Teuchos::rcp(new VectorSpaceThyra(model_->get_p_space(p_idx)));
	VectorSpace::space_ptr_t space_xD;
	//
	if(!no_model_x) {
		space_xD = Teuchos::rcp(new VectorSpaceThyra(model_space_x));
		if (p_idx >= 0)  {
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
	TEST_FOR_EXCEPT(!space_x_.get());

	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > model_space_f(model_->get_f_space());
	bool no_model_f = (model_space_f.get() == NULL);
	//
	if(!no_model_f)
		space_c_ = Teuchos::rcp(new VectorSpaceThyra(model_space_f));
	else
		space_c_ = Teuchos::null;
	//
	xinit_ = space_x_->create_member();  *xinit_ = 0.0;
	xl_    = space_x_->create_member();  *xl_    = -NLP::infinite_bound();
	xu_    = space_x_->create_member();  *xu_    = +NLP::infinite_bound();

	if(!no_model_f) {

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

		if(!no_model_x) {
			VectorSpace::vec_mut_ptr_t xinit_D = xinit_->sub_view(basis_sys_->var_dep());
			if(model_x0)  copy_from_model_x( model_x0, xinit_D.get() );
			else       copy_from_model_x( model_->get_x_init().get(), &*xinit_D );
			VectorSpace::vec_mut_ptr_t xl_D = xl_->sub_view(basis_sys_->var_dep());
			if(model_xL)  copy_from_model_x( model_xL, xl_D.get() );
			else       copy_from_model_x( model_->get_x_lower_bounds().get(), &*xl_D );
			VectorSpace::vec_mut_ptr_t xu_D = xu_->sub_view(basis_sys_->var_dep());
			if(model_xU)  copy_from_model_x( model_xU, xu_D.get() );
			else       copy_from_model_x( model_->get_x_upper_bounds().get(), &*xu_D );
		}

	}
  else {

		factory_Gc_ = Teuchos::null;

		basis_sys_ = Teuchos::null;

	}

  if(p_idx >= 0) {
    Range1D var_indep = ( basis_sys_.get() ? basis_sys_->var_indep() : Range1D() );
    VectorSpace::vec_mut_ptr_t xinit_I = xinit_->sub_view(var_indep);
    if(model_p0) copy_from_model_p( model_p0, &*xinit_I );
    else      copy_from_model_p( model_->get_p_init(p_idx).get(), &*xinit_I );
    VectorSpace::vec_mut_ptr_t xl_I = xl_->sub_view(var_indep);
    if(model_pL) copy_from_model_p( model_pL, &*xl_I );
    else      copy_from_model_p( model_->get_p_lower_bounds(p_idx).get(), &*xl_I );
    VectorSpace::vec_mut_ptr_t xu_I = xu_->sub_view(var_indep);
    if(model_pU) copy_from_model_p( model_pU, &*xu_I );
    else      copy_from_model_p( model_->get_p_upper_bounds(p_idx).get(), &*xu_I );
  }
  
	if(g_idx >= 0) {
		model_g_ = createMember(model_->get_g_space(g_idx));
	}

}

void NLPThyraModelEvaluator::assert_is_initialized() const
{
	TEST_FOR_EXCEPTION(
		!is_initialized(), NLP::UnInitialized
		,"NLPThyraModelEvaluator::assert_is_initialized() : Error, "
		"NLPThyraModelEvaluator::initialize() has not been called yet."
		);
}

void NLPThyraModelEvaluator::copy_from_model_x( const Thyra::VectorBase<value_type>* model_x, VectorMutable* x_D ) const
{
  if(!model_x) return;
	*x_D = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(model_x),false));
}

void NLPThyraModelEvaluator::copy_from_model_p( const Thyra::VectorBase<value_type> *model_p, VectorMutable* x_I ) const
{
  if(!model_p) return;
	*x_I = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(model_p),false));
}

void NLPThyraModelEvaluator::preprocessBaseInOutArgs(
  const Vector                                      &x
  ,bool                                             newx
  ,const ZeroOrderInfo                              *zero_order_info
  ,const ObjGradInfo                                *obj_grad_info
  ,const NLPFirstOrder::FirstOrderInfo              *first_order_info
  ,Thyra::ModelEvaluatorBase::InArgs<value_type>    *model_inArgs_inout
  ,Thyra::ModelEvaluatorBase::OutArgs<value_type>   *model_outArgs_inout
  ,MatrixOp*                                        *Gc_out
  ,VectorMutable*                                   *Gf_out
  ,value_type*                                      *f_out
  ,VectorMutable*                                   *c_out
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeMultiVector<value_type> DerivMV;
  typedef MEB::Derivative<value_type> Deriv;
  //
	if(newx) model_g_updated_ = model_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;
  //
  // Set the input arguments
  //
  MEB::InArgs<value_type> &model_inArgs = *model_inArgs_inout;
	if( basis_sys_.get() ) {
		const Range1D
			var_dep   = basis_sys_->var_dep(),
			var_indep = basis_sys_->var_indep();
		RefCountPtr<const Vector> xD = x.sub_view(var_dep), xI;
		if(p_idx_>=0) xI = x.sub_view(var_indep);
    model_inArgs.set_x(dyn_cast<const VectorMutableThyra>(*xD).thyra_vec());
		if(p_idx_ >= 0)
      model_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(*xI).thyra_vec());
  }
	else { // no dependent vars
    TEST_FOR_EXCEPT(p_idx_<0);
    model_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(x).thyra_vec());
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
  MEB::OutArgs<value_type> &model_outArgs = *model_outArgs_inout;
  if( f && (g_idx_>=0) && !f_updated_ ) {
    model_outArgs.set_g(g_idx_,model_g_); // ToDo: Make more general!
  }
  if( c && !c_updated_ ) {
    Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > thyra_c;
    get_thyra_vector(*space_c_,c,&thyra_c);
    model_outArgs.set_f(thyra_c);
  }
  if( Gf && !Gf_updated_ ) {
    if(g_idx_>=0) {
      if(p_idx_>=0) {
        const Range1D
          var_dep   = basis_sys_->var_dep(),
          var_indep = basis_sys_->var_indep();
        model_outArgs.set_DgDx(
          g_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::VectorBase<value_type> >(
              dyn_cast<VectorMutableThyra>(*Gf->sub_view(var_dep)).thyra_vec()
              )
            ,MEB::DERIV_TRANS_MV_BY_ROW
            )
          );
        model_outArgs.set_DgDp(
          g_idx_,p_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::VectorBase<value_type> >(
              dyn_cast<VectorMutableThyra>(*Gf->sub_view(var_indep)).thyra_vec()
              )
            ,MEB::DERIV_TRANS_MV_BY_ROW
            )
          );
      }
    }
  }
  if(Gc_out) *Gc_out = Gc;
  if(Gf_out) *Gf_out = Gf;
  if(f_out)  *f_out  = f;
  if(c_out)  *c_out  = c;
}

void NLPThyraModelEvaluator::postprocessBaseOutArgs(
  Thyra::ModelEvaluatorBase::OutArgs<value_type>        *model_outArgs_inout
  ,VectorMutable                                        *Gf
  ,value_type                                           *f
  ,VectorMutable                                        *c
  ) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  MEB::OutArgs<value_type> &model_outArgs = *model_outArgs_inout;
  if( f && !f_updated_ ) {
    if(g_idx_>=0) {
      *f = ::Thyra::get_ele(*model_g_,g_idx_);
    }
    else {
      *f = 0.0;
    }
    f_updated_ = true;
  }
  if( c && !c_updated_ ) {
    Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >
      thyra_c = model_outArgs.get_f();
    commit_thyra_vector(*space_c_,c,&thyra_c);
    model_outArgs.set_f(Teuchos::null);
    c_updated_ = true;
  }
  if( Gf && !Gf_updated_ ) {
    if(g_idx_>=0) {
      const Range1D
        var_dep   = basis_sys_->var_dep(),
        var_indep = basis_sys_->var_indep();
      Gf->sub_view(var_dep)->has_changed();
      Gf->sub_view(var_indep)->has_changed();
    }
    else {
      *Gf = 0.0;
    }
    Gf_updated_ = true;
  }
}

// private

void NLPThyraModelEvaluator::evalModel( 
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    *zero_order_info
  ,const ObjGradInfo      *obj_grad_info
  ) const
{
  typedef Thyra::ModelEvaluatorBase MEB;
  //
  // Set the input and output arguments
  //
  MEB::InArgs<value_type>  model_inArgs  = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  preprocessBaseInOutArgs(
    x,newx,zero_order_info,obj_grad_info,NULL
    ,&model_inArgs,&model_outArgs,NULL,&Gf,&f,&c
    );
  //
  // Evaluate the model
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the output arguments
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,c);
}

}	// end namespace NLPInterfacePack
