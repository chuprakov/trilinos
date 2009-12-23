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

#ifndef FEAPP_KL_EXPONENTIAL_FUNCTION_HPP
#define FEAPP_KL_EXPONENTIAL_FUNCTION_HPP

#include "FEApp_AbstractFunction.hpp"

#ifdef HAVE_SACADO_STOKHOS

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ParameterRegistration.hpp"
#include "Sacado_Traits.hpp"

#include "Stokhos_KL_ExponentialRandomField.hpp"

namespace FEApp {

  /*!
   * \brief A PDE material function based on a truncated KL expansion of
   * an exponential random field
   */
  template <typename EvalT>
  class KLExponentialFunction : 
    public FEApp::AbstractFunction<EvalT>,
    public Sacado::ParameterAccessor<EvalT, FEApp::EvaluationTraits> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractFunction<EvalT>::ScalarT ScalarT;
  
    //! Constructor
    KLExponentialFunction(Teuchos::ParameterList& params,
			  const Teuchos::RCP<ParamLib>& paramLib);

    //! Destructor
    virtual ~KLExponentialFunction() {};

    //! Evaluate function
    virtual void
    evaluate(const std::vector<double>& quad_points,
	     std::vector<ScalarT>& value) const;

    //! Accessor to the parameter for use by the Parameter Library
    ScalarT& getValue(const std::string &n);

  private:

    //! Private to prohibit copying
    KLExponentialFunction(const KLExponentialFunction&);

    //! Private to prohibit copying
    KLExponentialFunction& operator=(const KLExponentialFunction&);

  protected:

    //! Stokhos random field object
    Stokhos::KL::ExponentialRandomField<double> rf;

    //! Random variables
    Teuchos::Array<ScalarT> random_variables;

  };
}

// Include implementation
#include "FEApp_KLExponentialFunctionImpl.hpp"

#endif // HAVE_SACADO_STOKHOS

#endif // FEAPP_KL_EXPONENTIAL_FUNCTION_HPP
