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
#include "AbstractLinAlgPack_MatrixOpNonsingThyra.hpp"
#include "AbstractLinAlgPack_VectorSpaceBlocked.hpp"
#include "AbstractLinAlgPack_VectorAuxiliaryOps.hpp"
#include "AbstractLinAlgPack_BasisSystem.hpp"
#include "Thyra_ExplicitVectorView.hpp"
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
	initializeBase(model,p_idx,g_idx,model_xL,model_xU,model_x0,model_pL,model_pU,model_p0);
  TEST_FOR_EXCEPT(true);
}
	
// Overridden public members from NLP

void NLPDirectThyraModelEvaluator::initialize(bool test_setup)
{
  TEST_FOR_EXCEPT(true);
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
  TEST_FOR_EXCEPT(true);
}

void NLPDirectThyraModelEvaluator::calc_semi_newton_step(
  const Vector    &x
  ,VectorMutable  *c
  ,bool           recalc_c
  ,VectorMutable  *py
  ) const
{
  TEST_FOR_EXCEPT(true);
}

}	// end namespace NLPInterfacePack
