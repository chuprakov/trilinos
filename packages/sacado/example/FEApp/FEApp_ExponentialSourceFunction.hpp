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

#ifndef FEAPP_EXPONENTIALSOURCEFUNCTION_HPP
#define FEAPP_EXPONENTIALSOURCEFUNCTION_HPP

#include "FEApp_AbstractSourceFunction.hpp"

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ParameterRegistration.hpp"
#include "Sacado_Traits.hpp"


namespace FEApp {

  /*!
   * \brief An exponential PDE source function
   */
  template <typename EvalT>
  class ExponentialSourceFunction : 
    public FEApp::AbstractSourceFunction<EvalT>,
    public Sacado::ParameterAccessor<EvalT, FEApp::EvaluationTraits> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractSourceFunction<EvalT>::ScalarT ScalarT;
  
    //! Default constructor
    ExponentialSourceFunction(
	       const ScalarT& factor,
	       const Teuchos::RCP<ParamLib>& paramLib) : 
      alpha(factor) 
    {
      // Add nonlinear factor to parameter library
      std::string name = "Exponential Source Function Nonlinear Factor";
      new Sacado::ParameterRegistration<EvalT, FEApp::EvaluationTraits>(name, this, paramLib);
    };

    // !Destructor
    virtual ~ExponentialSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
             std::vector<ScalarT>& value) const {
      for (unsigned int i=0; i<solution.size(); i++) {
        value[i] = std::exp(alpha)*std::exp(solution[i]);
        //value[i] = -1.0;
      }
    };

    //! Parameter access method for use by the Parameter Library
    ScalarT& getValue(const std::string &n) {
      // For multiple parameters in the same file, can distinguish by string
      // if (n == "Exponential Source Function Nonlinear Factor") 
      return alpha;
    };

  private:

    //! Private to prohibit copying
    ExponentialSourceFunction(const ExponentialSourceFunction&);

    //! Private to prohibit copying
    ExponentialSourceFunction& operator=(const ExponentialSourceFunction&);

  protected:
  
    //! Factor
    ScalarT alpha;

  };
}

#endif // FEAPP_CUBICSOURCEFUNCTION_HPP
