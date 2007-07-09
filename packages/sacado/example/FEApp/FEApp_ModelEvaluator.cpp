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
  sacado_param_vec = Teuchos::rcp(new Sacado::ScalarParameterVector);
  if (param_names != Teuchos::null)
    app->getParamLib()->fillVector(*param_names, *sacado_param_vec);

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
  return 
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy, *(app->getJacobianGraph())));
}

EpetraExt::ModelEvaluator::InArgs
FEApp::ModelEvaluator::createInArgs() const
{
  InArgsSetup inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(IN_ARG_x,true);
  inArgs.set_Np(1); // 1 parameter vector
  if (app->isTransient()) {
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
  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
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
  Teuchos::RCP<Epetra_Vector> f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();

  // Cast W to a CrsMatrix, throw an exception if this fails
  Teuchos::RCP<Epetra_CrsMatrix> W_out_crs = 
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(W_out, true);
  
  //
  // Compute the functions
  //
  if(f_out != Teuchos::null && W_out != Teuchos::null) {
    app->computeGlobalJacobian(alpha, beta, x_dot, x, *sacado_param_vec,
			       *f_out, *W_out_crs);
  }
  else if(f_out != Teuchos::null ) {
    app->computeGlobalResidual(x_dot, x, *sacado_param_vec, *f_out);
  }
  else if(W_out != Teuchos::null ) {
    Epetra_Vector f(x.Map());
    app->computeGlobalJacobian(alpha, beta, x_dot, x, *sacado_param_vec, 
			       f, *W_out_crs);
  }
}
