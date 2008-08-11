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
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  template <typename ScalarT> class ConstantFunctionParameter;

  /*!
   * \brief A constant PDE function
   */
  template <typename ScalarT>
  class ConstantFunction : 
    public FEApp::AbstractFunction<ScalarT> {
  public:
  
    //! Default constructor
    ConstantFunction(
	       const ScalarT& value,
	       const Teuchos::RCP<Sacado::ScalarParameterLibrary>& paramLib) : 
      val(value) 
    {
      // Add val to parameter library
      std::string name = "Constant Function Value";
      if (!paramLib->isParameter(name))
	paramLib->addParameterFamily(name, true, false);
      if (!paramLib->template isParameterForType<ScalarT>(name)) {
	Teuchos::RCP< ConstantFunctionParameter<ScalarT> > tmp = 
	  Teuchos::rcp(new ConstantFunctionParameter<ScalarT>(Teuchos::rcp(this,false)));
	paramLib->template addEntry<ScalarT>(name, tmp);
      }
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
    }

    //! Set value
    void setValue(const ScalarT& value, bool mark_constant) { 
      val = value; 
      if (mark_constant) Sacado::MarkConstant<ScalarT>::eval(val); 
    }

    //! Get value
    const ScalarT& getValue() const { return val; }

  private:

    //! Private to prohibit copying
    ConstantFunction(const ConstantFunction&);

    //! Private to prohibit copying
    ConstantFunction& operator=(const ConstantFunction&);

  protected:
  
    //! Value
    ScalarT val;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * the value in the constant function
   */
  template <typename ScalarT>
  class ConstantFunctionParameter : 
    public Sacado::ScalarParameterEntry<ScalarT> {

  public:

    //! Constructor
    ConstantFunctionParameter(
		  const Teuchos::RCP< ConstantFunction<ScalarT> >& s) : 
      func(s) {}

    //! Destructor
    virtual ~ConstantFunctionParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      setValueAsConstant(ScalarT(value)); }

    //! Set parameter this object represents to \em value
    virtual void setValueAsConstant(const ScalarT& value) { 
      func->setValue(value, true); }
    
    //! Set parameter this object represents to \em value
    virtual void setValueAsIndependent(const ScalarT& value) { 
      func->setValue(value, false); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return func->getValue(); }
    
  protected:  
    
    //! Pointer to source function
    Teuchos::RCP< ConstantFunction<ScalarT> > func;

  };

}

#endif // FEAPP_CONSTANTFUNCTION_HPP
