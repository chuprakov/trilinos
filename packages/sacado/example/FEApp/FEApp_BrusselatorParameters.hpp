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

#ifndef FEAPP_BRUSSELATORPARAMETERS_HPP
#define FEAPP_BRUSSELATORPARAMETERS_HPP

#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"
#include "FEApp_TemplateTypes.hpp"

namespace FEApp {

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * alpha in the brusselator problem
   */
  template <typename EvalT>
  class BrusselatorAlphaParameter : 
    public Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> {

  public:

    //! Scalar type
    typedef typename Sacado::ScalarParameterEntry<EvalT,EvaluationTraits>::ScalarT ScalarT;

    //! Constructor
    BrusselatorAlphaParameter(const ScalarT& value) : val(value) {}

    //! Destructor
    virtual ~BrusselatorAlphaParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      val = value; Sacado::MarkConstant<ScalarT>::eval(val); }
    
    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      val = value; }

    //! Get real parameter value
    virtual double getRealValue() const {
      return Sacado::ScalarValue<ScalarT>::eval(val); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return val; }
    
  protected:  
    
    //! Value
    ScalarT val;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * beta in the brusselator problem
   */
  template <typename EvalT>
  class BrusselatorBetaParameter : 
    public Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> {

  public:

    //! Scalar type
    typedef typename Sacado::ScalarParameterEntry<EvalT,EvaluationTraits>::ScalarT ScalarT;

    //! Constructor
    BrusselatorBetaParameter(const ScalarT& value) : val(value) {}

    //! Destructor
    virtual ~BrusselatorBetaParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      val = value; Sacado::MarkConstant<ScalarT>::eval(val); }
    
    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      val = value; }

    //! Get real parameter value
    virtual double getRealValue() const {
      return Sacado::ScalarValue<ScalarT>::eval(val); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return val; }
    
  protected:  
    
    //! Value
    ScalarT val;

  };

}

#endif // FEAPP_BRUSSELATORPARAMETRS_HPP
