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
  const Teuchos::RCP< Teuchos::Array<std::string> >& free_param_names,
  const Teuchos::RCP< Teuchos::Array<std::string> >& sg_param_names) 
  : app(app_),
    supports_sg(false)
{
  // Compute number of parameter vectors
  int num_param_vecs = 1;
#if SG_ACTIVE
  if (sg_param_names != Teuchos::null) {
    supports_sg = true;
    num_param_vecs = 2;
  }
#endif

  param_names.resize(num_param_vecs);
  sacado_param_vec.resize(num_param_vecs);
  epetra_param_map.resize(num_param_vecs);
  epetra_param_vec.resize(num_param_vecs);

  // Set parameter names
  param_names[0] = free_param_names;
  if (num_param_vecs == 2)
    param_names[1] = sg_param_names;

  // Initialize each parameter vector
  const Epetra_Comm& comm = app->getMap()->Comm();
  for (int i=0; i<num_param_vecs; i++) {

    // Initialize Sacado parameter vector
    sacado_param_vec[i] = Teuchos::rcp(new ParamVec);
    if (param_names[i] != Teuchos::null)
      app->getParamLib()->fillVector<FEApp::ResidualType>(*(param_names[i]), 
							  *(sacado_param_vec[i]));

    // Create Epetra map for parameter vector
    epetra_param_map[i] = 
      Teuchos::rcp(new Epetra_LocalMap(sacado_param_vec[i]->size(), 0, comm));

    // Create Epetra vector for parameters
    epetra_param_vec[i] = 
      Teuchos::rcp(new Epetra_Vector(*(epetra_param_map[i])));
  
    // Set parameters
    for (unsigned int j=0; j<sacado_param_vec[i]->size(); j++)
      (*(epetra_param_vec[i]))[j] = (*(sacado_param_vec[i]))[j].baseValue;

  }

  // Create storage for SG parameter values
#if SG_ACTIVE
  if (supports_sg)
    p_sg_vals.resize(sg_param_names->size());
#endif
						      
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
  TEST_FOR_EXCEPTION(l >= static_cast<int>(epetra_param_map.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return epetra_param_map[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::ModelEvaluator::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return param_names[l];
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_x_init() const
{
  return app->getInitialSolution();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return epetra_param_vec[l];
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
  inArgs.set_Np(param_names.size());
  if (supports_sg) {
    inArgs.setSupports(IN_ARG_x_sg,true);
    inArgs.set_Np_sg(1); // 1 SG parameter vector
  }
  else
    inArgs.set_Np_sg(0);
  if (app->isTransient()) {
    inArgs.setSupports(IN_ARG_t,true);
    inArgs.setSupports(IN_ARG_x_dot,true);
    inArgs.setSupports(IN_ARG_alpha,true);
    inArgs.setSupports(IN_ARG_beta,true);
    if (supports_sg)
      inArgs.setSupports(IN_ARG_x_dot_sg,true);
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
  if (supports_sg) {
    outArgs.setSupports(OUT_ARG_f_sg,true);
    outArgs.setSupports(OUT_ARG_W_sg,true);
  }

  return outArgs;
}

void 
FEApp::ModelEvaluator::evalModel(const InArgs& inArgs, 
				 const OutArgs& outArgs) const
{
  //
  // Get the input arguments
  //
  Teuchos::RCP<const Epetra_Vector> x = inArgs.get_x();
  Teuchos::RCP<const Epetra_Vector> x_dot;
  double alpha = 0.0;
  double beta = 1.0;
  if (app->isTransient()) {
    x_dot = inArgs.get_x_dot();
    if (x_dot != Teuchos::null) {
      alpha = inArgs.get_alpha();
      beta = inArgs.get_beta();
    }
  }
  for (int i=0; i<inArgs.Np(); i++) {
    Teuchos::RCP<const Epetra_Vector> p = inArgs.get_p(i);
    if (p != Teuchos::null) {
      for (unsigned int j=0; j<sacado_param_vec[i]->size(); j++)
	(*(sacado_param_vec[i]))[j].baseValue = (*p)[j];
    }
  }

  //
  // Get the output arguments
  //
  EpetraExt::ModelEvaluator::Evaluation<Epetra_Vector> f_out = outArgs.get_f();
  Teuchos::RCP<Epetra_Operator> W_out = outArgs.get_W();
  
  //
  // Compute the functions
  //
  bool f_computed = false;

  // W matrix
  if (W_out != Teuchos::null) {
    if (f_out.getType() == EVAL_TYPE_EXACT ||
        f_out.getType() == EVAL_TYPE_APPROX_DERIV)
      app->computeGlobalJacobian(alpha, beta, x_dot.get(), *x, 
				 sacado_param_vec,
                                 f_out.get(), *W_out);
    else
      app->computeGlobalPreconditioner(alpha, beta, x_dot.get(), *x, 
                                       sacado_param_vec, 
				       f_out.get(), *W_out);
    f_computed = true;
  }
  
  // df/dp
  if (outArgs.Np() > 0) {
    for (int i=0; i<outArgs.Np(); i++) {
      Teuchos::RCP<Epetra_MultiVector> dfdp_out = 
	outArgs.get_DfDp(i).getMultiVector();
      if (dfdp_out != Teuchos::null) {
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DfDp(i).getDerivativeMultiVector().getParamIndexes();
	unsigned int n_params = p_indexes.size();
	Teuchos::Array<std::string> p_names(n_params);
	for (unsigned int j=0; j<n_params; j++)
	  p_names[j] = (*(param_names[i]))[p_indexes[j]];
	ParamVec p_vec;
	app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, p_vec);
	for (unsigned int j=0; j<p_vec.size(); j++)
	  p_vec[j].baseValue = (*(sacado_param_vec[i]))[p_indexes[j]].baseValue;
  
	app->computeGlobalTangent(0.0, 0.0, false, x_dot.get(), *x, 
				  sacado_param_vec, &p_vec,
				  NULL, NULL, f_out.get(), NULL, 
				  dfdp_out.get());

	f_computed = true;
      }
    }
  }

  // f
  if(f_out != Teuchos::null && !f_computed) {
    app->computeGlobalResidual(x_dot.get(), *x, sacado_param_vec, 
			       *f_out);
  }

#if SG_ACTIVE
  if (supports_sg) {
    InArgs::sg_const_vector_t x_sg = 
      inArgs.get_x_sg();
    if (x_sg != Teuchos::null) {
      InArgs::sg_const_vector_t x_dot_sg;
      if (app->isTransient())
	x_dot_sg = inArgs.get_x_dot_sg();
      InArgs::sg_const_vector_t epetra_p_sg = 
	inArgs.get_p_sg(0);
      Teuchos::Array<SGType> *p_sg_ptr = NULL;
      if (epetra_p_sg != Teuchos::null) {
	for (unsigned int i=0; i<p_sg_vals.size(); i++) {
	  int num_sg_blocks = epetra_p_sg->size();
	  p_sg_vals[i].resize(num_sg_blocks);
	  p_sg_vals[i].copyForWrite();
	  for (int j=0; j<num_sg_blocks; j++) {
	    p_sg_vals[i].fastAccessCoeff(j) = (*epetra_p_sg)[j][i];
	  }
	}
	p_sg_ptr = &p_sg_vals;
      }
      OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
      OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
      if (W_sg != Teuchos::null)
	app->computeGlobalSGJacobian(alpha, beta, x_dot_sg.get(), *x_sg, 
				     sacado_param_vec[0].get(), 
				     sacado_param_vec[1].get(), p_sg_ptr,
				     f_sg.get(), *W_sg);
      else if (f_sg != Teuchos::null)
	app->computeGlobalSGResidual(x_dot_sg.get(), *x_sg, 
				     sacado_param_vec[0].get(), 
				     sacado_param_vec[1].get(), p_sg_ptr,
				     *f_sg);
    }
  }
#endif
}
