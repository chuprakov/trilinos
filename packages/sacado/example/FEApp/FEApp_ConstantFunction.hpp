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

#ifndef FEAPP_CONSTANTFUNCTION_HPP
#define FEAPP_CONSTANTFUNCTION_HPP

#include "FEApp_AbstractFunction.hpp"

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ParameterRegistration.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  /*!
   * \brief A constant PDE function
   */
  template <typename EvalT>
  class ConstantFunction : 
    public FEApp::AbstractFunction<EvalT>,
    public Sacado::ParameterAccessor<EvalT, FEApp::EvaluationTraits> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractFunction<EvalT>::ScalarT ScalarT;
  
    //! Default constructor
    ConstantFunction(
	       const ScalarT& value,
	       const Teuchos::RCP<ParamLib>& paramLib) : 
      val(value) 
    {
      // Add val to parameter library
      std::string name = "Constant Function Value";
      // Register this name and *this* class with the Parameter Library
      new Sacado::ParameterRegistration<EvalT, FEApp::EvaluationTraits>(name, this, paramLib);
    };

    //! Destructor
    virtual ~ConstantFunction() {};

    //! Evaluate function
    virtual void
    evaluate(const std::vector<double>& quad_points,
	     std::vector<ScalarT>& value) const {
      for (unsigned int i=0; i<quad_points.size(); i++) {
        value[i] = val;
      }
    };

    //! Accessor to the parameter for use by the Parameter Library
    ScalarT& getValue(const std::string &n) {
      if (n == "Constant Function Value") return val;
      else {
        TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                       std::endl <<
                       "Error! Logic error in setting paramter " << n
                       << " in FEApp_ConstantFunction.hpp" << std::endl)
        return val;
      }
    };

  private:

    //! Private to prohibit copying
    ConstantFunction(const ConstantFunction&);

    //! Private to prohibit copying
    ConstantFunction& operator=(const ConstantFunction&);

  protected:
  
    //! Value
    ScalarT val;

  };
}

#endif // FEAPP_CONSTANTFUNCTION_HPP
