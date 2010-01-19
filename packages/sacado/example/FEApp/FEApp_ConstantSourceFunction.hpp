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

#ifndef FEAPP_CONSTANTSOURCEFUNCTION_HPP
#define FEAPP_CONSTANTSOURCEFUNCTION_HPP

#include "FEApp_AbstractSourceFunction.hpp"

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  template <typename EvalT> class ConstantSourceFunctionParameter;

  /*!
   * \brief A quadratic PDE source function
   */
  template <typename EvalT>
  class ConstantSourceFunction : 
    public FEApp::AbstractSourceFunction<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractSourceFunction<EvalT>::ScalarT ScalarT;
  
    //! Default constructor
    ConstantSourceFunction(
	       const ScalarT& value,
	       const Teuchos::RCP<ParamLib>& paramLib) : 
      v(value) 
    {
      // Add nonlinear factor to parameter library
      std::string name = "Constant Source Function Value";
      if (!paramLib->isParameter(name))
        paramLib->addParameterFamily(name, true, false);
      if (!paramLib->template isParameterForType<EvalT>(name)) {
        Teuchos::RCP< ConstantSourceFunctionParameter<EvalT> > tmp = 
          Teuchos::rcp(new ConstantSourceFunctionParameter<EvalT>(Teuchos::rcp(this,false)));
        paramLib->template addEntry<EvalT>(name, tmp);
      }
    };

    //! Destructor
    virtual ~ConstantSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
             std::vector<ScalarT>& value) const {
      for (unsigned int i=0; i<solution.size(); i++)
        value[i] = v;
    }

    //! Set nonlinear factor
    void setValue(const ScalarT& val, bool mark_constant) { 
      v = val; 
      if (mark_constant) Sacado::MarkConstant<ScalarT>::eval(v); 
    }

    //! Get nonlinear factor
    const ScalarT& getValue() const { return v; }

  private:

    //! Private to prohibit copying
    ConstantSourceFunction(const ConstantSourceFunction&);

    //! Private to prohibit copying
    ConstantSourceFunction& operator=(const ConstantSourceFunction&);

  protected:
  
    //! Value
    ScalarT v;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * the value of the constant source function
   */
  template <typename EvalT>
  class ConstantSourceFunctionParameter : 
    public Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> {

  public:

    //! Scalar type
    typedef typename Sacado::ScalarParameterEntry<EvalT,EvaluationTraits>::ScalarT ScalarT;

    //! Constructor
    ConstantSourceFunctionParameter(
		  const Teuchos::RCP< ConstantSourceFunction<EvalT> >& s) : 
      srcFunc(s) {}

    //! Destructor
    virtual ~ConstantSourceFunctionParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      srcFunc->setValue(value, true); }
    
    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      srcFunc->setValue(value, false); }

     //! Get real parameter value
    virtual double getRealValue() const { 
      return Sacado::ScalarValue<ScalarT>::eval(srcFunc->getValue()); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return srcFunc->getValue(); }
    
  protected:  
    
    //! Pointer to source function
    Teuchos::RCP< ConstantSourceFunction<EvalT> > srcFunc;

  };

}

#endif // FEAPP_CONSTANTSOURCEFUNCTION_HPP
