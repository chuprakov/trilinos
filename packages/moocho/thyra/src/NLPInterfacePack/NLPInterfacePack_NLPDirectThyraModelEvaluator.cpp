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

#include "NLPInterfacePack_NLPDirectThyraModelEvaluator.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "AbstractLinAlgPack_VectorOut.hpp"
#include "AbstractLinAlgPack_ThyraAccessors.hpp"
#include "AbstractLinAlgPack_VectorSpaceThyra.hpp"
#include "AbstractLinAlgPack_VectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MultiVectorMutableThyra.hpp"
#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "AbstractLinAlgPack_LinAlgOpPack.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {

NLPDirectThyraModelEvaluator::NLPDirectThyraModelEvaluator()
{}

NLPDirectThyraModelEvaluator::NLPDirectThyraModelEvaluator(
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
	initialize(model,p_idx,g_idx,model_xL,model_xU,model_x0,model_pL,model_pU,model_p0);
}

void NLPDirectThyraModelEvaluator::initialize(
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
  typedef Thyra::ModelEvaluatorBase MEB;
	initializeBase(model,p_idx,g_idx,model_xL,model_xU,model_x0,model_pL,model_pU,model_p0);
  Thyra::ModelEvaluatorBase::OutArgs<double> model_outArgs = model->createOutArgs();
  MEB::DerivativeProperties model_W_properties = model_outArgs.get_W_properties();
  if( p_idx >= 0 ) {
    TEST_FOR_EXCEPTION(
      !model_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).supports(MEB::DERIV_MV_BY_COL),std::invalid_argument
      ,"Error, model must support computing DfDp("<<p_idx<<") as a column-oriented multi-vector!"
      );
  }
}

// Overridden public members from NLP

void NLPDirectThyraModelEvaluator::initialize(bool test_setup)
{
	if(initialized_) {
		NLPDirect::initialize(test_setup);
		return;
	}
  NLPThyraModelEvaluator::initialize(test_setup);
  NLPDirect::initialize(test_setup);
}

void NLPDirectThyraModelEvaluator::unset_quantities()
{
  NLPDirect::unset_quantities();
}

// Overridden public members from NLPDirect

Range1D NLPDirectThyraModelEvaluator::var_dep() const
{
  return basis_sys_->var_dep();
}

Range1D NLPDirectThyraModelEvaluator::var_indep() const
{
  return basis_sys_->var_indep();
}

const NLPDirect::mat_fcty_ptr_t
NLPDirectThyraModelEvaluator::factory_D() const
{
  return basis_sys_->factory_D();
}

const NLPDirect::mat_sym_nonsing_fcty_ptr_t
NLPDirectThyraModelEvaluator::factory_S() const
{
  return basis_sys_->factory_S();
}

void NLPDirectThyraModelEvaluator::calc_point(
  const Vector     &x
  ,value_type      *f
  ,VectorMutable   *c
  ,bool            recalc_c
  ,VectorMutable   *Gf
  ,VectorMutable   *py
  ,VectorMutable   *rGf
  ,MatrixOp        *GcU
  ,MatrixOp        *D
  ,MatrixOp        *Uz
  ) const
{
  using Teuchos::dyn_cast;
  using Teuchos::RefCountPtr;
  using Teuchos::rcp;
  using Teuchos::rcp_const_cast;
  using Teuchos::rcp_dynamic_cast;
  using AbstractLinAlgPack::VectorSpaceThyra;
  using AbstractLinAlgPack::VectorMutableThyra;
  using AbstractLinAlgPack::MultiVectorMutableThyra;
  using AbstractLinAlgPack::MatrixOpThyra;
  using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef Thyra::ModelEvaluatorBase MEB;
  typedef MEB::DerivativeMultiVector<value_type> DerivMV;
  typedef MEB::Derivative<value_type> Deriv;
  //
  // Validate input
  //
  TEST_FOR_EXCEPT(GcU!=NULL);  // Can't handle these yet!
  TEST_FOR_EXCEPT(Uz!=NULL);
  //
  // Set the input and output arguments
  //
  // ToDo: Disallow computation Gf and instead just compute
  // the operators for this!
  //
  MEB::InArgs<value_type>  model_inArgs  = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  NLPObjGrad::ObjGradInfo obj_grad_info;
  obj_grad_info.Gf = Gf;
  obj_grad_info.f = f;
  if(recalc_c) obj_grad_info.c = c;
  preprocessBaseInOutArgs(
    x,true,NULL,&obj_grad_info,NULL
    ,&model_inArgs,&model_outArgs,NULL,NULL,NULL,NULL
    );
  if( py || rGf || D ) {
    if(thyra_C_.get()==NULL)
      thyra_C_ = model_->create_W();
    model_outArgs.set_W(thyra_C_);
  }
  if( rGf || D ) {
    if(thyra_N_.get()==NULL)
      thyra_N_ = model_->create_DfDp_mv(p_idx_,MEB::DERIV_MV_BY_COL).getMultiVector();
    model_outArgs.set_DfDp(p_idx_,DerivMV(thyra_N_,MEB::DERIV_MV_BY_COL));
  }
  //
  // Evaluate the functions
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the evaluation
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,recalc_c?c:NULL);
  // Setup solve components
  const VectorSpaceThyra                                       *space_c;
  const VectorSpaceThyra                                       *space_xD;
  Teuchos::RefCountPtr<const Thyra::VectorBase<value_type> >   thyra_c;
  Teuchos::RefCountPtr<Thyra::VectorBase<value_type> >         thyra_py;
  RefCountPtr<MatrixOp>                                        D_used = rcp(D,false);
  RefCountPtr<Thyra::MultiVectorBase<value_type> >             thyra_D;
  if( py ) {
    space_c  = &dyn_cast<const VectorSpaceThyra>(c->space()),
    space_xD = &dyn_cast<const VectorSpaceThyra>(py->space());
    get_thyra_vector(*space_c,*c,&thyra_c);
    get_thyra_vector(*space_xD,py,&thyra_py);
  }
  if( D || rGf ) {
    if(!D) D_used = this->factory_D()->create();
    thyra_D =
      rcp_const_cast<Thyra::MultiVectorBase<value_type> >(
        rcp_dynamic_cast<const Thyra::MultiVectorBase<value_type> >(
          dyn_cast<MultiVectorMutableThyra>(*D_used).thyra_multi_vec()
          )
        );
  }
  // Perform solve
  if( ( D || rGf ) && py ) {
    // Solve for [py,D] all at once!
    const int nind = thyra_N_->domain()->dim();
    RefCountPtr<Thyra::MultiVectorBase<value_type> > 
      thyra_cN = Thyra::createMembers(thyra_N_->range(),nind+1);
    Thyra::assign(&*thyra_cN->col(0),*thyra_c);
    Thyra::assign(&*thyra_cN->subView(Teuchos::Range1D(1,nind)),*thyra_N_);
    RefCountPtr<Thyra::MultiVectorBase<value_type> > 
      thyra_pyD = Thyra::createMembers(thyra_D->range(),nind+1);
    Thyra::assign(&*thyra_pyD,0.0);
    Thyra::solve(*thyra_C_,Thyra::NOTRANS,*thyra_cN,&*thyra_pyD);
    Thyra::scale(-1.0,&*thyra_pyD);
    Thyra::assign(&*thyra_py,*thyra_pyD->col(0));
    Thyra::assign(&*thyra_D,*thyra_pyD->subView(Teuchos::Range1D(1,nind)));
  }
  else {
    // Solve for py or D
    if( py ) {
      // py = -inv(C)*c
      Thyra::assign(&*thyra_py,0.0);
      Thyra::solve(*thyra_C_,Thyra::NOTRANS,*thyra_c,&*thyra_py);
      Thyra::Vt_S(&*thyra_py,-1.0);
    }
    if( D || rGf ) {
      // D = -inv(C)*N
      Thyra::assign(&*thyra_D,0.0);
      Thyra::solve(*thyra_C_,Thyra::NOTRANS,*thyra_N_,&*thyra_D);
      Thyra::scale(-1.0,&*thyra_D);
      // ToDo: Just compute the operators allocated with Gf and not Gf directly!
    }
  }
  if(thyra_py.get()) {
    free_thyra_vector(*space_c,*c,&thyra_c);
    commit_thyra_vector(*space_xD,py,&thyra_py);
  }
  // Compute reduced gradient
  if(rGf) {
    // rGf = D' * Gf_xD + Gf_xI
    const Range1D
      var_dep   = basis_sys_->var_dep(),
      var_indep = basis_sys_->var_indep();
    LinAlgOpPack::V_MtV( rGf, *D_used, BLAS_Cpp::trans, *Gf->sub_view(var_dep) );
    LinAlgOpPack::Vp_V( rGf, *Gf->sub_view(var_indep) );
  }
  // * ToDo: Add specialized algorithm for computing D using an inexact Jacobian
  // * ToDo: Add in logic for inexact solves
}

void NLPDirectThyraModelEvaluator::calc_semi_newton_step(
  const Vector    &x
  ,VectorMutable  *c
  ,bool           recalc_c
  ,VectorMutable  *py
  ) const
{
  if(recalc_c) {
    //
    // Recompute c
    //
    TEST_FOR_EXCEPT(true);
  }
  // Compute py = - inv(C)*c
  TEST_FOR_EXCEPT(true);
}

}	// end namespace NLPInterfacePack
