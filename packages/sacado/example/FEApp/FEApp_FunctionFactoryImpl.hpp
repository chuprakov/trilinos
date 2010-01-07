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
#include "Teuchos_TestForException.hpp"
#include "FEApp_ConstantFunction.hpp"
#include "FEApp_KLExponentialFunction.hpp"

template <typename EvalT>
FEApp::FunctionFactory<EvalT>::FunctionFactory(
	    const Teuchos::RCP<Teuchos::ParameterList>& funcParams_,
	    const Teuchos::RCP<ParamLib>& paramLib_) :
  funcParams(funcParams_), paramLib(paramLib_)
{
}

template <typename EvalT>
Teuchos::RCP< FEApp::AbstractFunction<EvalT> >
FEApp::FunctionFactory<EvalT>::create()
{
  Teuchos::RCP< FEApp::AbstractFunction<EvalT> > strategy;

  std::string& method = funcParams->get("Name", "Constant");
  if (method == "Constant") {
    double value = funcParams->get("Constant Value", 1.0);
    strategy = 
      Teuchos::rcp(new FEApp::ConstantFunction<EvalT>(value, paramLib));
  }
#ifdef HAVE_SACADO_STOKHOS
  else if (method == "KL Exponential Random Field")
    strategy = 
      Teuchos::rcp(new FEApp::KLExponentialFunction<EvalT>(*funcParams, 
							   paramLib));
#endif
  else {
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                       std::endl << 
                       "Error!  Unknown function " << method << 
                       "!" << std::endl << "Supplied parameter list is " << 
                       std::endl << *funcParams);
  }

  return strategy;
}
