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

// Debugging only
//#include <iostream>
//#include <typeinfo>
//#include "AbstractLinAlgPack_VectorOut.hpp"
//#include "TSFCoreTestingTools.hpp"

/** ToDo:

11/28/2005: Add some flags to determine if we are going to be using a direct
or adjoint sensitivities.  If we are using direct sensitivities we want
np_dfdp to be a multi-vector and if we are using adjoint sensitivities we want
np_dfdp to be a general linear operator.

11/28/2005: Change the interface to only allow a single objective function and
a single set of parameter vectors.  That will greatly simplify things and that
is all the MOOCHO/TSFCore interface did anyway.  We can always build decorator
subclasses of Thyra::ModelEvaluator to build different types of transformed
models.

*/

namespace NLPInterfacePack {

NLPThyraModelEvaluator::NLPThyraModelEvaluator()
	:initialized_(false),obj_scale_(1.0)
{}

NLPThyraModelEvaluator::NLPThyraModelEvaluator(
  const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np   
  ,const int                                                      p_idx 
  ,const int                                                      g_idx 
  ,const Thyra::VectorBase<value_type>                            *np_xL
  ,const Thyra::VectorBase<value_type>                            *np_xU
  ,const Thyra::VectorBase<value_type>                            *np_x0
  ,const Thyra::VectorBase<value_type>                            *np_pL
  ,const Thyra::VectorBase<value_type>                            *np_pU
  ,const Thyra::VectorBase<value_type>                            *np_p0
	)
	:initialized_(false),obj_scale_(1.0)
{
	initialize(np,p_idx,g_idx,np_xL,np_xU,np_x0,np_pL,np_pU,np_p0);
}

void NLPThyraModelEvaluator::initialize(
  const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >  &np
  ,const int                                                      p_idx
  ,const int                                                      g_idx
  ,const Thyra::VectorBase<value_type>                            *np_xL
  ,const Thyra::VectorBase<value_type>                            *np_xU
  ,const Thyra::VectorBase<value_type>                            *np_x0
  ,const Thyra::VectorBase<value_type>                            *np_pL
  ,const Thyra::VectorBase<value_type>                            *np_pU
  ,const Thyra::VectorBase<value_type>                            *np_p0
	)
{

	using Teuchos::dyn_cast;
	using AbstractLinAlgPack::VectorSpaceThyra;
	using AbstractLinAlgPack::VectorMutableThyra;
	using AbstractLinAlgPack::MatrixOpNonsingThyra;
  typedef ::Thyra::ModelEvaluatorBase MEB;
	
	initialized_ = false;
	np_g_updated_ = np_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;

	const char msg_err[] = "NLPThyraModelEvaluator::initialize(...): Errror!";
  Thyra::ModelEvaluatorBase::OutArgs<double> np_outArgs = np->createOutArgs();
	TEST_FOR_EXCEPTION( np.get() == NULL, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( p_idx >= 0 && ( p_idx == 0 || p_idx > np_outArgs.Np() ), std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( g_idx >= 0 && ( g_idx == 0 || g_idx > np_outArgs.Ng() ), std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( !np_outArgs.supports(MEB::OUT_ARG_W), std::invalid_argument, msg_err );
  MEB::DerivativeProperties np_W_properties = np_outArgs.get_W_properties();
	TEST_FOR_EXCEPTION( np_W_properties.supportsAdjoint==false, std::invalid_argument, msg_err );
	TEST_FOR_EXCEPTION( np_W_properties.rank==MEB::DERIV_RANK_DEFICIENT, std::invalid_argument, msg_err );
  if(p_idx >= 0 ) {
    TEST_FOR_EXCEPTION( np_outArgs.supports(MEB::OUT_ARG_DfDp,p_idx).none(), std::invalid_argument, msg_err );
    if(g_idx >= 0) {
      TEST_FOR_EXCEPTION( np_outArgs.supports(MEB::OUT_ARG_DgDp,g_idx,p_idx).none(), std::invalid_argument, msg_err );
    }
  }
  if(g_idx >= 0) {
    TEST_FOR_EXCEPTION( np_outArgs.supports(MEB::OUT_ARG_DgDx,g_idx).none(), std::invalid_argument, msg_err );
  }
	//if(!np->isInitialized()) np->initialize();
	//
	np_ = np;
  p_idx_ = p_idx;
  g_idx_ = g_idx;
	//
	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > np_space_x(np_->get_x_space());
	bool no_np_x = (np_space_x.get() == NULL);
	//
	VectorSpace::space_ptr_t space_xI;
  if(p_idx >= 0)
    space_xI = Teuchos::rcp(new VectorSpaceThyra(np_->get_p_space(p_idx)));
	VectorSpace::space_ptr_t space_xD;
	//
	if(!no_np_x) {
		space_xD = Teuchos::rcp(new VectorSpaceThyra(np_space_x));
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

	Teuchos::RefCountPtr<const Thyra::VectorSpaceBase<value_type> > np_space_f(np_->get_f_space());
	bool no_np_f = (np_space_f.get() == NULL);
	//
	if(!no_np_f)
		space_c_ = Teuchos::rcp(new VectorSpaceThyra(np_space_f));
	else
		space_c_ = Teuchos::null;
	//
	xinit_ = space_x_->create_member();  *xinit_ = 0.0;
	xl_    = space_x_->create_member();  *xl_    = -NLP::infinite_bound();
	xu_    = space_x_->create_member();  *xu_    = +NLP::infinite_bound();

	if(!no_np_f) {

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

		if(!no_np_x) {
			VectorSpace::vec_mut_ptr_t xinit_D = xinit_->sub_view(basis_sys_->var_dep());
			if(np_x0)  copy_from_np_x( np_x0, xinit_D.get() );
			else       copy_from_np_x( np_->get_x_init().get(), &*xinit_D );
			VectorSpace::vec_mut_ptr_t xl_D = xl_->sub_view(basis_sys_->var_dep());
			if(np_xL)  copy_from_np_x( np_xL, xl_D.get() );
			else       copy_from_np_x( np_->get_x_lower_bounds().get(), &*xl_D );
			VectorSpace::vec_mut_ptr_t xu_D = xu_->sub_view(basis_sys_->var_dep());
			if(np_xU)  copy_from_np_x( np_xU, xu_D.get() );
			else       copy_from_np_x( np_->get_x_upper_bounds().get(), &*xu_D );
		}

	}
  else {

		factory_Gc_ = Teuchos::null;

		basis_sys_ = Teuchos::null;

	}

  if(p_idx >= 0) {
    Range1D var_indep = ( basis_sys_.get() ? basis_sys_->var_indep() : Range1D() );
    VectorSpace::vec_mut_ptr_t xinit_I = xinit_->sub_view(var_indep);
    if(np_p0) copy_from_np_p( np_p0, &*xinit_I );
    else      copy_from_np_p( np_->get_p_init(p_idx).get(), &*xinit_I );
    VectorSpace::vec_mut_ptr_t xl_I = xl_->sub_view(var_indep);
    if(np_pL) copy_from_np_p( np_pL, &*xl_I );
    else      copy_from_np_p( np_->get_p_lower_bounds(p_idx).get(), &*xl_I );
    VectorSpace::vec_mut_ptr_t xu_I = xu_->sub_view(var_indep);
    if(np_pU) copy_from_np_p( np_pU, &*xu_I );
    else      copy_from_np_p( np_->get_p_upper_bounds(p_idx).get(), &*xu_I );
  }
  
	if(g_idx >= 0) {
		np_g_ = createMember(np_->get_g_space(g_idx));
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

// Overridden public members from NLPFirstOrder

void NLPThyraModelEvaluator::set_Gc(MatrixOp* Gc)
{
  NLPFirstOrder::set_Gc(Gc);
  Gc_updated_ = false;
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

void NLPThyraModelEvaluator::copy_from_np_x( const Thyra::VectorBase<value_type>* np_x, VectorMutable* x_D ) const
{
  if(!np_x) return;
	*x_D = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(np_x),false));
}

void NLPThyraModelEvaluator::copy_from_np_p( const Thyra::VectorBase<value_type> *np_p, VectorMutable* x_I ) const
{
  if(!np_p) return;
	*x_I = AbstractLinAlgPack::VectorMutableThyra(Teuchos::rcp(const_cast<Thyra::VectorBase<value_type>*>(np_p),false));
}

void NLPThyraModelEvaluator::evalModel( 
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
	if(newx) np_g_updated_ = np_Dg_updated_ = f_updated_ = c_updated_ = Gf_updated_ = Gc_updated_ = false;
  //
  // Set the input arguments
  //
  MEB::InArgs<value_type>
    np_inArgs = np_->createInArgs();
	if( basis_sys_.get() != NULL) {
		const Range1D
			var_dep   = basis_sys_->var_dep(),
			var_indep = basis_sys_->var_indep();
		RefCountPtr<const Vector> xD = x.sub_view(var_dep), xI;
		if(p_idx_>=0) xI = x.sub_view(var_indep);
    np_inArgs.set_x(dyn_cast<const VectorMutableThyra>(*xD).thyra_vec());
		if(p_idx_ >= 0)
      np_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(*xI).thyra_vec());
  }
	else { // no dependent vars
    TEST_FOR_EXCEPT(p_idx_ < 0);
    np_inArgs.set_p(p_idx_,dyn_cast<const VectorMutableThyra>(x).thyra_vec());
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
  MEB::OutArgs<value_type> np_outArgs = np_->createOutArgs();
  if( f && (g_idx_>=0) && !f_updated_ ) {
    np_outArgs.set_g(g_idx_,np_g_); // ToDo: Make more general!
  }
  if( c && !c_updated_ ) {
    Teuchos::RefCountPtr<Thyra::VectorBase<value_type> > thyra_c;
    get_thyra_vector(*space_c_,c,&thyra_c);
    np_outArgs.set_f(thyra_c);
  }
  if( Gf && !Gf_updated_ ) {
    if(g_idx_>=0) {
      if(p_idx_>=0) {
        const Range1D
          var_dep   = basis_sys_->var_dep(),
          var_indep = basis_sys_->var_indep();
        np_outArgs.set_DgDx(
          g_idx_
          ,DerivMV(
            rcp_const_cast<Thyra::VectorBase<value_type> >(
              dyn_cast<VectorMutableThyra>(*Gf->sub_view(var_dep)).thyra_vec()
              )
            ,MEB::DERIV_TRANS_MV_BY_ROW
            )
          );
        np_outArgs.set_DgDp(
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

  MatrixOpNonsing  *C_aggr;
  MatrixOp         *N_aggr;
  if( Gc && !Gc_updated_ ) {
		BasisSystemComposite::get_C_N( Gc, &C_aggr, &N_aggr ); // Will return NULLs if Gc is not initialized
    if(C_aggr) {
      np_outArgs.set_W(
        rcp_const_cast<Thyra::LinearOpWithSolveBase<value_type> >(
          dyn_cast<MatrixOpNonsingThyra>(*C_aggr).set_uninitialized()
          )
        );
      if(p_idx_ >= 0) {
        // ToDo: This is implemented for direct sensitivities, update for adjoints also!
        np_outArgs.set_DfDp(
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
      np_outArgs.set_W(np_->create_W());
      if(p_idx_>=0)
        np_outArgs.set_DfDp(p_idx_,np_->create_DfDp_mv(1,MEB::DERIV_MV_BY_COL));
    }
    if(np_inArgs.supports(MEB::IN_ARG_alpha)) np_inArgs.set_alpha(0.0);
    np_inArgs.set_beta(1.0);
  }
  //
  // Evaluate the model
  //
  np_->evalModel(np_inArgs,np_outArgs);
  //
  // Update outputs
  //
  if( f && !f_updated_ ) {
    if(g_idx_>=0) {
      *f = ::Thyra::get_ele(*np_g_,g_idx_);
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
    dyn_cast<MatrixOpNonsingThyra>(*C_aggr).initialize(np_outArgs.get_W(),BLAS_Cpp::no_trans);
    if(p_idx_>=0)
      dyn_cast<MatrixOpThyra>(*N_aggr).initialize(np_outArgs.get_DfDp(p_idx_).getDerivativeMultiVector().getMultiVector(),BLAS_Cpp::no_trans);
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
