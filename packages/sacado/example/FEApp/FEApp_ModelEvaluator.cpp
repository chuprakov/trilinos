// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "FEApp_ModelEvaluator.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"

FEApp::ModelEvaluator::ModelEvaluator(
  const Teuchos::RCP<FEApp::Application>& app_,
  const Teuchos::RCP< Teuchos::Array<std::string> >& free_param_names) 
  : app(app_),
    param_names(free_param_names)
{
  // Initialize Sacado parameter vector
  sacado_param_vec = Teuchos::rcp(new ParamVec);
  if (param_names != Teuchos::null)
    app->getParamLib()->fillVector<FEApp::ResidualType>(*param_names, 
                                                        *sacado_param_vec);

  // Create Epetra map for parameter vector
  const Epetra_Comm& comm = app->getMap()->Comm();
  epetra_param_map = Teuchos::rcp(new Epetra_LocalMap(sacado_param_vec->size(),
                                                      0, comm));

  // Create Epetra vector for parameters
  epetra_param_vec = Teuchos::rcp(new Epetra_Vector(*epetra_param_map));
  
  // Set parameters
  for (unsigned int i=0; i<sacado_param_vec->size(); i++)
    (*epetra_param_vec)[i] = (*sacado_param_vec)[i].baseValue;
						      
}

// Overridden from EpetraExt::ModelEvaluator

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_x_map() const
{
  return app->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_f_map() const
{
  return app->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_p_map(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_map() only " <<
                     " supports 1 parameter vector.  Supplied index l = " << 
                     l << std::endl);

  return epetra_param_map;
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::ModelEvaluator::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names() only " <<
                     " supports 1 parameter vector.  Supplied index l = " << 
                     l << std::endl);
  return param_names;
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_x_init() const
{
  return app->getInitialSolution();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_init() only " <<
                     " supports 1 parameter vector.  Supplied index l = " << 
                     l << std::endl);
  
  return epetra_param_vec;
}

Teuchos::RCP<Epetra_Operator>
FEApp::ModelEvaluator::create_W() const
{
  return app->createW();
}

EpetraExt::ModelEvaluator::InArgs
FEApp::ModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1); // 1 parameter vector
  if (app->isTransient()) {
    inArgs.setSupports(IN_ARG_t,true);
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
  }
  return inArgs;
}

EpetraExt::ModelEvaluator::OutArgs
FEApp::ModelEvaluator::createOutArgs() const
{
  OutArgsSetup outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.set_Np_Ng(1, 0); // 1 parameter vector
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.setSupports(OUT_ARG_DfDp, 0, DerivativeSupport(DERIV_MV_BY_COL));
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL ,true)
    );
  return outArgs;
}

void 
FEApp::ModelEvaluator::evalModel(const InArgs& inArgs, 
				 const OutArgs& outArgs) const
{
  //
  // Get the input arguments
  //
  const Epetra_Vector& x = *inArgs.get_x();
  const Epetra_Vector *x_dot = NULL;
  double alpha = 0.0;
  double beta = 1.0;
  if (app->isTransient()) {
    x_dot = inArgs.get_x_dot().get();
    if (x_dot) {
      alpha = inArgs.get_alpha();
      beta = inArgs.get_beta();
    }
  }
  const Epetra_Vector *p = inArgs.get_p(0).get();
  if (p != NULL) {
    for (unsigned int i=0; i<sacado_param_vec->size(); i++)
      (*sacado_param_vec)[i].baseValue = (*p)[i];
  }

  //
  // Get the output arguments
  //
  //Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  Teuchos::RCP<Epetra_MultiVector> dfdp_out;
  if (outArgs.Np() > 0)
    dfdp_out = outArgs.get_DfDp(0).getMultiVector();
  
  //
  // Compute the functions
  //
  if(W_out != Teuchos::null) {
    if (f_out.getType() == EVAL_TYPE_EXACT ||
        f_out.getType() == EVAL_TYPE_APPROX_DERIV)
      app->computeGlobalJacobian(alpha, beta, x_dot, x, sacado_param_vec.get(),
                                 f_out.get(), *W_out);
    else
      app->computeGlobalPreconditioner(alpha, beta, x_dot, x, 
                                       sacado_param_vec.get(), f_out.get(), 
                                       *W_out);
  }
  else if (dfdp_out != Teuchos::null) {
    Teuchos::Array<int> p_indexes = 
      outArgs.get_DfDp(0).getDerivativeMultiVector().getParamIndexes();
    unsigned int n_params = p_indexes.size();
    Teuchos::Array<std::string> p_names(n_params);
    for (unsigned int i=0; i<n_params; i++)
      p_names[i] = (*param_names)[p_indexes[i]];
    ParamVec p_vec;
    app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, p_vec);
    for (unsigned int i=0; i<p_vec.size(); i++)
      p_vec[i].baseValue = (*p)[p_indexes[i]];
  
    app->computeGlobalTangent(0.0, 0.0, false, x_dot, x, &p_vec,
                              NULL, NULL, f_out.get(), NULL, dfdp_out.get());
  }
  else if(f_out != Teuchos::null ) {
    app->computeGlobalResidual(x_dot, x, sacado_param_vec.get(), *f_out);
  }
}
