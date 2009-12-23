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

template <typename EvalT>
FEApp::KLExponentialFunction<EvalT>::
KLExponentialFunction(Teuchos::ParameterList& params,
		      const Teuchos::RCP<ParamLib>& paramLib) :
  rf(params),
  random_variables(rf.stochasticDimension())
{
  // Add random variables to parameter library
  int dim = rf.stochasticDimension();
  for (int i=0; i<dim; i++) {
    std::stringstream name;
    name << "KL Exponential Function Random Variable " << i;
    new Sacado::ParameterRegistration<EvalT, FEApp::EvaluationTraits>(name.str(), this, paramLib);
  }

  // Initialize random variables
  for (int i=0; i<dim; i++) {
    std::stringstream ss;
    ss << "KL Exponential Function Random Variable " << i;
    random_variables[i] = params.get(ss.str(), 0.0);
  }
}

template <typename EvalT>
void
FEApp::KLExponentialFunction<EvalT>::
evaluate(const std::vector<double>& quad_points,
	 std::vector<ScalarT>& value) const {
  Teuchos::Array<double> point(1);
  for (unsigned int i=0; i<quad_points.size(); i++) {
    point[0] = quad_points[i];
    value[i] = rf.evaluate(point, random_variables);
  }
}

template <typename EvalT>
typename FEApp::KLExponentialFunction<EvalT>::ScalarT&
FEApp::KLExponentialFunction<EvalT>::
getValue(const std::string &n) 
{
  int dim = rf.stochasticDimension();
  for (int i=0; i<dim; i++) {
    std::stringstream name;
    name << "KL Exponential Function Random Variable " << i;
    if (n == name.str()) 
      return random_variables[i];
  }
  TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
		     std::endl <<
		     "Error! Logic error in setting paramter " << n
		     << " in FEApp_KLExponentialFunction.hpp" << std::endl);
  return random_variables[0];
}
