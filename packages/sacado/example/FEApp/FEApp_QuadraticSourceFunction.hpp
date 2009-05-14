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

#ifndef FEAPP_QUADRATICSOURCEFUNCTION_HPP
#define FEAPP_QUADRATICSOURCEFUNCTION_HPP

#include "FEApp_AbstractSourceFunction.hpp"

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  template <typename EvalT> class QuadraticNonlinearFactorParameter;

  /*!
   * \brief A quadratic PDE source function
   */
  template <typename EvalT>
  class QuadraticSourceFunction : 
    public FEApp::AbstractSourceFunction<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractSourceFunction<EvalT>::ScalarT ScalarT;
  
    //! Default constructor
    QuadraticSourceFunction(
	       const ScalarT& factor,
	       const Teuchos::RCP<ParamLib>& paramLib) : 
      alpha(factor) 
    {
      // Add nonlinear factor to parameter library
      std::string name = "Quadratic Source Function Nonlinear Factor";
      if (!paramLib->isParameter(name))
        paramLib->addParameterFamily(name, true, false);
      if (!paramLib->template isParameterForType<EvalT>(name)) {
        Teuchos::RCP< QuadraticNonlinearFactorParameter<EvalT> > tmp = 
          Teuchos::rcp(new QuadraticNonlinearFactorParameter<EvalT>(Teuchos::rcp(this,false)));
        paramLib->template addEntry<EvalT>(name, tmp);
      }
    };

    //! Destructor
    virtual ~QuadraticSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
             std::vector<ScalarT>& value) const {
      for (unsigned int i=0; i<solution.size(); i++)
        value[i] = alpha*solution[i]*solution[i];
    }

    //! Set nonlinear factor
    void setFactor(const ScalarT& val, bool mark_constant) { 
      alpha = val; 
      if (mark_constant) Sacado::MarkConstant<ScalarT>::eval(alpha); 
    }

    //! Get nonlinear factor
    const ScalarT& getFactor() const { return alpha; }

  private:

    //! Private to prohibit copying
    QuadraticSourceFunction(const QuadraticSourceFunction&);

    //! Private to prohibit copying
    QuadraticSourceFunction& operator=(const QuadraticSourceFunction&);

  protected:
  
    //! Factor
    ScalarT alpha;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * the nonlinear factor in the quadratic source function
   */
  template <typename EvalT>
  class QuadraticNonlinearFactorParameter : 
    public Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> {

  public:

    //! Scalar type
    typedef typename Sacado::ScalarParameterEntry<EvalT,EvaluationTraits>::ScalarT ScalarT;

    //! Constructor
    QuadraticNonlinearFactorParameter(
		  const Teuchos::RCP< QuadraticSourceFunction<EvalT> >& s) : 
      srcFunc(s) {}

    //! Destructor
    virtual ~QuadraticNonlinearFactorParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      srcFunc->setFactor(value, true); }
    
    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      srcFunc->setFactor(value, false); }

     //! Get real parameter value
    virtual double getRealValue() const { 
      return Sacado::ScalarValue<ScalarT>::eval(srcFunc->getFactor()); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return srcFunc->getFactor(); }
    
  protected:  
    
    //! Pointer to source function
    Teuchos::RCP< QuadraticSourceFunction<EvalT> > srcFunc;

  };

}

#endif // FEAPP_QUADRATICSOURCEFUNCTION_HPP
