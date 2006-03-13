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

#include <assert.h>

#include <algorithm>

#include "NLPInterfacePack_NLPFirstOrderThyraModelEvaluator.hpp"
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
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace NLPInterfacePack {

NLPFirstOrderThyraModelEvaluator::NLPFirstOrderThyraModelEvaluator()
{}

NLPFirstOrderThyraModelEvaluator::NLPFirstOrderThyraModelEvaluator(
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

void NLPFirstOrderThyraModelEvaluator::initialize(
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
	initializeBase(model,p_idx,g_idx,model_xL,model_xU,model_x0,model_pL,model_pU,model_p0);
}
	
// Overridden public members from NLP

void NLPFirstOrderThyraModelEvaluator::initialize(bool test_setup)
{
	if(initialized_) {
		NLPFirstOrder::initialize(test_setup);
		return;
	}
  NLPThyraModelEvaluator::initialize(test_setup);
  NLPFirstOrder::initialize(test_setup);
}

void NLPFirstOrderThyraModelEvaluator::unset_quantities()
{
  NLPFirstOrder::unset_quantities();
}

// Overridden public members from NLPFirstOrder

void NLPFirstOrderThyraModelEvaluator::set_Gc(MatrixOp* Gc)
{
  NLPFirstOrder::set_Gc(Gc);
  Gc_updated_ = false;
}

const NLPFirstOrder::mat_fcty_ptr_t
NLPFirstOrderThyraModelEvaluator::factory_Gc() const
{
	return factory_Gc_;
}

const NLPFirstOrder::basis_sys_ptr_t
NLPFirstOrderThyraModelEvaluator::basis_sys() const
{
	return basis_sys_;
}

// Overridden protected members from NLPFirstOrder

void NLPFirstOrderThyraModelEvaluator::imp_calc_Gc(const Vector& x, bool newx, const FirstOrderInfo& first_order_info) const
{
  evalModel(x,newx,NULL,NULL,&first_order_info);
}

// private

void NLPFirstOrderThyraModelEvaluator::evalModel( 
  const Vector            &x
  ,bool                   newx
  ,const ZeroOrderInfo    *zero_order_info
  ,const ObjGradInfo      *obj_grad_info
  ,const FirstOrderInfo   *first_order_info
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
  // Set the input and output arguments
  //
  MEB::InArgs<value_type>  model_inArgs  = model_->createInArgs();
  MEB::OutArgs<value_type> model_outArgs = model_->createOutArgs();
  MatrixOp            *Gc = NULL;
  VectorMutable       *Gf = NULL;
  value_type          *f  = NULL;
  VectorMutable       *c  = NULL;
  preprocessBaseInOutArgs(
    x,newx,zero_order_info,obj_grad_info,first_order_info
    ,&model_inArgs,&model_outArgs,&Gc,&Gf,&f,&c
    );
  //
  MatrixOpNonsing  *C_aggr;
  MatrixOp         *N_aggr;
  if( Gc && !Gc_updated_ ) {
		BasisSystemComposite::get_C_N( Gc, &C_aggr, &N_aggr ); // Will return NULLs if Gc is not initialized
    if(C_aggr) {
      model_outArgs.set_W(
        rcp_const_cast<Thyra::LinearOpWithSolveBase<value_type> >(
          dyn_cast<MatrixOpNonsingThyra>(*C_aggr).set_uninitialized()
          )
        );
      if(p_idx_ >= 0) {
        // ToDo: This is implemented for direct sensitivities, change this for adjoint sensitivities
        model_outArgs.set_DfDp(
          p_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::MultiVectorBase<value_type> >(
              rcp_dynamic_cast<const Thyra::MultiVectorBase<value_type> >(
                dyn_cast<MatrixOpThyra>(*N_aggr).set_uninitialized()
                )
              )
            ,MEB::DERIV_MV_BY_COL
            )
          );
      }
    }
    else {
      model_outArgs.set_W(model_->create_W());
      if(p_idx_>=0)
        model_outArgs.set_DfDp(p_idx_,model_->create_DfDp_mv(p_idx_,MEB::DERIV_MV_BY_COL));
    }
    if(model_inArgs.supports(MEB::IN_ARG_alpha)) model_inArgs.set_alpha(0.0);
  }
  //
  // Evaluate the model
  //
  model_->evalModel(model_inArgs,model_outArgs);
  //
  // Postprocess the output arguments
  //
  postprocessBaseOutArgs(&model_outArgs,Gf,f,c);
  //
  if( Gc && !Gc_updated_ ) {
    RefCountPtr<MatrixOpNonsing> C_ptr;
    RefCountPtr<MatrixOp>        N_ptr;
    if(!C_aggr) {
      C_ptr  = Teuchos::rcp(new MatrixOpNonsingThyra());
      C_aggr = &*C_ptr;
      if(p_idx_>=0) {
        N_ptr  = Teuchos::rcp(new MatrixOpThyra());
        N_aggr = &*N_ptr;
      }
    }
    dyn_cast<MatrixOpNonsingThyra>(*C_aggr).initialize(model_outArgs.get_W(),BLAS_Cpp::no_trans);
    // ToDo: This is implemented for direct sensitivities, change this for adjoint sensitivities
    if(p_idx_>=0)
      dyn_cast<MatrixOpThyra>(*N_aggr).initialize(model_outArgs.get_DfDp(p_idx_).getDerivativeMultiVector().getMultiVector(),BLAS_Cpp::no_trans);
    if( C_ptr.get() ) {
      BasisSystemComposite::initialize_Gc(
        this->space_x(), basis_sys_->var_dep(), basis_sys_->var_indep()
        ,this->space_c()
        ,C_ptr, N_ptr
        ,Gc
        );
    }
    Gc_updated_ = true;
  }
}

}	// end namespace NLPInterfacePack
