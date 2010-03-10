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
  const Teuchos::RCP< Teuchos::Array<std::string> >& sg_param_names
#if SG_ACTIVE
  ,const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>& initial_x_sg_,
  const Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>& initial_p_sg_
#endif
) 
  : app(app_),
    supports_p(false),
    supports_g(false),
    supports_sg(false),
#if SG_ACTIVE
    initial_x_sg(initial_x_sg_),
    initial_p_sg(initial_p_sg_),
#endif
    eval_W_with_f(false)
{
  // Compute number of parameter vectors
  int num_param_vecs = 0;
  if (free_param_names != Teuchos::null) {
    supports_p = true;
    num_param_vecs = 1;
  }
#if SG_ACTIVE
  if (sg_param_names != Teuchos::null) {
    supports_sg = true;
    num_param_vecs = 2;
  }
#endif

  if (num_param_vecs > 0) {
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

  supports_g = (app->getResponseMap() != Teuchos::null);
						      
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
  TEST_FOR_EXCEPTION(supports_p == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_map():  " <<
                     "No parameters have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l >= static_cast<int>(epetra_param_map.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return epetra_param_map[l];
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_p_sg_map(int l) const
{
  TEST_FOR_EXCEPTION(supports_sg == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_map():  " <<
                     "SG is not enabled.");
  TEST_FOR_EXCEPTION(l != 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_map():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return epetra_param_map[1];
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_g_map(int l) const
{
  TEST_FOR_EXCEPTION(supports_g == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_map():  " <<
                     "No response functions have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_map() only " <<
                     " supports 1 response vector.  Supplied index l = " << 
                     l << std::endl);

  return app->getResponseMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::ModelEvaluator::get_g_sg_map(int l) const
{
  TEST_FOR_EXCEPTION(supports_sg == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_sg_map():  " <<
                     "SG is not enabled.");
  TEST_FOR_EXCEPTION(l != 0, Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_g_sg_map() only " <<
                     " supports 1 response vector.  Supplied index l = " << 
                     l << std::endl);

  return app->getResponseMap();
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::ModelEvaluator::get_p_names(int l) const
{
  TEST_FOR_EXCEPTION(supports_p == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names():  " <<
                     "No parameters have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return param_names[l];
}

Teuchos::RCP<const Teuchos::Array<std::string> >
FEApp::ModelEvaluator::get_p_sg_names(int l) const
{
  TEST_FOR_EXCEPTION(supports_sg == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_names():  " <<
                     "SG is not enabled.");
  TEST_FOR_EXCEPTION(l != 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_names():  " <<
                     "Invalid parameter index l = " << l << std::endl);

  return param_names[1];
}

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_x_init() const
{
  return app->getInitialSolution();
}

#if SG_ACTIVE
Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>
FEApp::ModelEvaluator::get_x_sg_init() const
{
  TEST_FOR_EXCEPTION(supports_sg == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_x_sg_init():  " <<
                     "SG is not enabled.");
  TEST_FOR_EXCEPTION(initial_x_sg == Teuchos::null, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_x_sg_init():  " <<
                     "initial_x_sg is NULL!");
  return initial_x_sg;
}
#endif

Teuchos::RCP<const Epetra_Vector>
FEApp::ModelEvaluator::get_p_init(int l) const
{
  TEST_FOR_EXCEPTION(supports_p == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_init():  " <<
                     "No parameters have been supplied.  " <<
                     "Supplied index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(l >= static_cast<int>(param_names.size()) || l < 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  
  return epetra_param_vec[l];
}

#if SG_ACTIVE
Teuchos::RCP<const Stokhos::EpetraVectorOrthogPoly>
FEApp::ModelEvaluator::get_p_sg_init(int l) const
{
  TEST_FOR_EXCEPTION(supports_sg == false, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_init():  " <<
                     "SG is not enabled.");
  TEST_FOR_EXCEPTION(l != 0, 
		     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_init():  " <<
                     "Invalid parameter index l = " << l << std::endl);
  TEST_FOR_EXCEPTION(initial_p_sg == Teuchos::null, 
                     Teuchos::Exceptions::InvalidParameter,
                     std::endl << 
                     "Error!  FEApp::ModelEvaluator::get_p_sg_init():  " <<
                     "initial_x_sg is NULL!");
  return initial_p_sg;
}
#endif

Teuchos::RCP<Epetra_Operator>
FEApp::ModelEvaluator::create_W() const
{
  my_W = app->createW();
  return my_W;
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
    inArgs.setSupports(IN_ARG_sg_basis,true);
    inArgs.setSupports(IN_ARG_sg_quadrature,true);
    inArgs.setSupports(IN_ARG_sg_expansion,true);
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

  if (supports_p && supports_g) {
    outArgs.set_Np_Ng(param_names.size(), 1);
    for (int i=0; i<param_names.size(); i++)
      outArgs.setSupports(OUT_ARG_DgDp, 0, i, 
			  DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
  }
  else if (supports_p)
    outArgs.set_Np_Ng(1, 0);
  else if (supports_g)
    outArgs.set_Np_Ng(0, 1);
  else
    outArgs.set_Np_Ng(0, 0);

  if (supports_g) {
    outArgs.setSupports(OUT_ARG_DgDx, 0, 
			DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    if (app->isTransient())
      outArgs.setSupports(OUT_ARG_DgDx_dot, 0, 
                          DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
  }
  if (supports_p)
    for (int i=0; i<param_names.size(); i++)
      outArgs.setSupports(OUT_ARG_DfDp, i, 
			  DerivativeSupport(DERIV_TRANS_MV_BY_ROW));

  outArgs.setSupports(OUT_ARG_f,true);
  outArgs.setSupports(OUT_ARG_W,true);
  outArgs.set_W_properties(
    DerivativeProperties(
      DERIV_LINEARITY_UNKNOWN ,DERIV_RANK_FULL ,true)
    );

  if (supports_sg) {
    outArgs.setSupports(OUT_ARG_f_sg,true);
    outArgs.setSupports(OUT_ARG_W_sg,true);

    if (supports_p && supports_g) {
      outArgs.set_Np_Ng_sg(param_names.size(), 1);
      for (int i=0; i<param_names.size(); i++)
	outArgs.setSupports(OUT_ARG_DgDp_sg, 0, i, 
			    DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    }
    else if (supports_p)
      outArgs.set_Np_Ng_sg(1, 0);
    else if (supports_g)
      outArgs.set_Np_Ng_sg(0, 1);
    else
      outArgs.set_Np_Ng_sg(0, 0);

    if (supports_g) {
      outArgs.setSupports(OUT_ARG_DgDx_sg, 0, 
			  DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
      if (app->isTransient())
	outArgs.setSupports(OUT_ARG_DgDx_dot_sg, 0, 
			    DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
    }
    if (supports_p)
      for (int i=0; i<param_names.size(); i++)
	outArgs.setSupports(OUT_ARG_DfDp_sg, i, 
			    DerivativeSupport(DERIV_TRANS_MV_BY_ROW));
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
  if (f_out != Teuchos::null && eval_W_with_f) {
    app->computeGlobalJacobian(alpha, beta, x_dot.get(), *x, sacado_param_vec,
                                 f_out.get(), *my_W);
    f_computed = true;
  }
  else if (W_out != Teuchos::null && !eval_W_with_f) {
    if (f_out.getType() == EVAL_TYPE_EXACT ||
        f_out.getType() == EVAL_TYPE_APPROX_DERIV)
      app->computeGlobalJacobian(alpha, beta, x_dot.get(), *x, sacado_param_vec,
                                 f_out.get(), *W_out);
    else
      app->computeGlobalPreconditioner(alpha, beta, x_dot.get(), *x, 
                                       sacado_param_vec, 
				       f_out.get(), *W_out);
    f_computed = true;
  }
  
  // df/dp
  if (supports_p) {
    for (int i=0; i<outArgs.Np(); i++) {
      Teuchos::RCP<Epetra_MultiVector> dfdp_out = 
	outArgs.get_DfDp(i).getMultiVector();
      if (dfdp_out != Teuchos::null) {
	Teuchos::Array<int> p_indexes = 
	  outArgs.get_DfDp(i).getDerivativeMultiVector().getParamIndexes();
	unsigned int n_params = p_indexes.size();
	Teuchos::RCP<ParamVec> p_vec;
	if (n_params > 0) {
	  Teuchos::Array<std::string> p_names(n_params);
	  for (unsigned int j=0; j<n_params; j++)
	    p_names[j] = (*(param_names[i]))[p_indexes[j]];
	  p_vec = Teuchos::rcp(new ParamVec);
	  app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, *p_vec);
	  for (unsigned int j=0; j<p_vec->size(); j++)
	    (*p_vec)[j].baseValue = 
	      (*(sacado_param_vec[i]))[p_indexes[j]].baseValue;
	}
	else
	  p_vec = Teuchos::rcp(sacado_param_vec[0].get(), false);
  
	app->computeGlobalTangent(0.0, 0.0, false, x_dot.get(), *x, 
				  sacado_param_vec, p_vec.get(),
				  NULL, NULL, f_out.get(), NULL, 
				  dfdp_out.get());

	f_computed = true;
      }
    }
  }

  // f
  if(f_out != Teuchos::null && !f_computed) {
    app->computeGlobalResidual(x_dot.get(), *x, sacado_param_vec, *f_out);
  }

  // Response functions
  if (outArgs.Ng() > 0 && supports_g) {
    Teuchos::RCP<Epetra_Vector> g_out = outArgs.get_g(0);
    Teuchos::RCP<Epetra_MultiVector> dgdx_out = 
      outArgs.get_DgDx(0).getMultiVector();
    Teuchos::RCP<Epetra_MultiVector> dgdxdot_out;
    if (app->isTransient())
      dgdxdot_out = outArgs.get_DgDx_dot(0).getMultiVector();
    
    Teuchos::Array< Teuchos::RCP<ParamVec> > p_vec(outArgs.Np());
    Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dgdp_out(outArgs.Np());
    bool have_dgdp = false;
    for (int i=0; i<outArgs.Np(); i++) {
      dgdp_out[i] = outArgs.get_DgDp(0,i).getMultiVector();
      if (dgdp_out[i] != Teuchos::null)
	have_dgdp = true;
      Teuchos::Array<int> p_indexes = 
	outArgs.get_DgDp(0,i).getDerivativeMultiVector().getParamIndexes();
      unsigned int n_params = p_indexes.size();
      if (n_params > 0) {
	Teuchos::Array<std::string> p_names(n_params);
	for (unsigned int j=0; j<n_params; j++)
	  p_names[i] = (*(param_names[i]))[p_indexes[j]];
	p_vec[i] = Teuchos::rcp(new ParamVec);
	app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, 
							    *(p_vec[i]));
	for (unsigned int j=0; j<p_vec[i]->size(); j++)
	  (*(p_vec[i]))[j].baseValue = (*(sacado_param_vec[i]))[p_indexes[j]].baseValue;
      }
      else
	p_vec[i] = sacado_param_vec[i];
    }

    if (have_dgdp ||dgdx_out != Teuchos::null || dgdxdot_out != Teuchos::null) {
      app->evaluateResponseGradients(x_dot.get(), *x, sacado_param_vec, p_vec, 
                                     g_out.get(), dgdx_out.get(), 
                                     dgdxdot_out.get(), dgdp_out);
    }
    else if (g_out != Teuchos::null)
      app->evaluateResponses(x_dot.get(), *x, sacado_param_vec, *g_out);
  }

#if SG_ACTIVE
  if (supports_sg) {
    InArgs::sg_const_vector_t x_sg = 
      inArgs.get_x_sg();
    if (x_sg != Teuchos::null) {
      app->init_sg(inArgs.get_sg_basis(), inArgs.get_sg_quadrature(), 
		   inArgs.get_sg_expansion());
      InArgs::sg_const_vector_t x_dot_sg;
      if (app->isTransient())
	x_dot_sg = inArgs.get_x_dot_sg();
      InArgs::sg_const_vector_t epetra_p_sg = 
	inArgs.get_p_sg(0);
      Teuchos::Array<SGType> *p_sg_ptr = NULL;
      if (epetra_p_sg != Teuchos::null) {
	for (int i=0; i<p_sg_vals.size(); i++) {
	  int num_sg_blocks = epetra_p_sg->size();
	  p_sg_vals[i].reset(app->getStochasticExpansion(), num_sg_blocks);
	  p_sg_vals[i].copyForWrite();
	  for (int j=0; j<num_sg_blocks; j++) {
	    p_sg_vals[i].fastAccessCoeff(j) = (*epetra_p_sg)[j][i];
	  }
	}
	p_sg_ptr = &p_sg_vals;
      }

      OutArgs::sg_vector_t f_sg = outArgs.get_f_sg();
      OutArgs::sg_operator_t W_sg = outArgs.get_W_sg();
      bool f_sg_computed = false;
      
      // W_sg
      if (W_sg != Teuchos::null) {
	app->computeGlobalSGJacobian(alpha, beta, x_dot_sg.get(), *x_sg, 
				     sacado_param_vec[0].get(), 
				     sacado_param_vec[1].get(), p_sg_ptr,
				     f_sg.get(), *W_sg);
	f_sg_computed = true;
      }

      // df/dp_sg
      if (supports_p) {
	for (int i=0; i<outArgs.Np_sg(); i++) {
	  Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dfdp_sg 
	    = outArgs.get_DfDp_sg(i).getMultiVector();
	  if (dfdp_sg != Teuchos::null) {
	    Teuchos::Array<int> p_indexes = 
	      outArgs.get_DfDp_sg(i).getDerivativeMultiVector().getParamIndexes();
	    unsigned int n_params = p_indexes.size();
	    Teuchos::RCP<ParamVec> p_vec;
	    if (n_params > 0) {
	      Teuchos::Array<std::string> p_names(n_params);
	      for (unsigned int j=0; j<n_params; j++)
		p_names[j] = (*(param_names[i]))[p_indexes[j]];
	      p_vec = Teuchos::rcp(new ParamVec);
	      app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, *p_vec);
	      for (unsigned int j=0; j<p_vec->size(); j++)
		(*p_vec)[j].baseValue = 
		  (*(sacado_param_vec[i]))[p_indexes[j]].baseValue;
	    }
	    else
	      p_vec = Teuchos::rcp(sacado_param_vec[0].get(), false);
	    
	    app->computeGlobalSGTangent(0.0, 0.0, false, x_dot_sg.get(), *x_sg, 
					sacado_param_vec[0].get(), p_vec.get(), 
					sacado_param_vec[1].get(), p_sg_ptr,
					NULL, NULL, f_sg.get(), NULL, 
					dfdp_sg.get());

	    f_sg_computed = true;
	  }
	}
      }

      if (f_sg != Teuchos::null && !f_sg_computed)
	app->computeGlobalSGResidual(x_dot_sg.get(), *x_sg, 
				     sacado_param_vec[0].get(), 
				     sacado_param_vec[1].get(), p_sg_ptr,
				     *f_sg);

      // Response functions
      if (outArgs.Ng_sg() > 0 && supports_g) {
	Teuchos::RCP< Stokhos::EpetraVectorOrthogPoly > g_sg 
	  = outArgs.get_g_sg(0);
	Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdx_sg 
	  = outArgs.get_DgDx_sg(0).getMultiVector();
	Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > dgdxdot_sg;
	if (app->isTransient())
	  dgdxdot_sg = outArgs.get_DgDx_dot_sg(0).getMultiVector();
    
	Teuchos::Array< Teuchos::RCP<ParamVec> > p_vec(outArgs.Np());
	Teuchos::Array< Teuchos::RCP< Stokhos::EpetraMultiVectorOrthogPoly > > dgdp_sg(outArgs.Np());
	bool have_dgdp = false;
	for (int i=0; i<outArgs.Np(); i++) {
	  dgdp_sg[i] = outArgs.get_DgDp_sg(0,i).getMultiVector();
	  if (dgdp_sg[i] != Teuchos::null)
	    have_dgdp = true;
	  Teuchos::Array<int> p_indexes = 
	    outArgs.get_DgDp_sg(0,i).getDerivativeMultiVector().getParamIndexes();
	  unsigned int n_params = p_indexes.size();
	  if (n_params > 0) {
	    Teuchos::Array<std::string> p_names(n_params);
	    for (unsigned int j=0; j<n_params; j++)
	      p_names[i] = (*(param_names[i]))[p_indexes[j]];
	    p_vec[i] = Teuchos::rcp(new ParamVec);
	    app->getParamLib()->fillVector<FEApp::ResidualType>(p_names, 
								*(p_vec[i]));
	    for (unsigned int j=0; j<p_vec[i]->size(); j++)
	      (*(p_vec[i]))[j].baseValue = (*(sacado_param_vec[i]))[p_indexes[j]].baseValue;
	  }
	  else
	    p_vec[i] = sacado_param_vec[i];
	}

	if (have_dgdp ||dgdx_sg != Teuchos::null || 
	    dgdxdot_sg != Teuchos::null) {
	  app->evaluateSGResponseGradients(x_dot_sg.get(), *x_sg, 
					   sacado_param_vec, p_vec, p_sg_ptr,
					   g_sg.get(), dgdx_sg.get(), 
					   dgdxdot_sg.get(), dgdp_sg);
	}
	else if (g_sg != Teuchos::null)
	  app->evaluateSGResponses(x_dot_sg.get(), *x_sg, 
				   sacado_param_vec, p_sg_ptr, *g_sg);
      }
    }
  }
#endif
}
