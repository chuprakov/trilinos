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
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  template <typename ScalarT> class ExponentialNonlinearFactorParameter;

  /*!
   * \brief A cubic PDE source function
   */
  template <typename ScalarT>
  class ExponentialSourceFunction : 
    public FEApp::AbstractSourceFunction<ScalarT> {
  public:
  
    //! Default constructor
    ExponentialSourceFunction(
	       const ScalarT& factor,
	       const Teuchos::RCP<Sacado::ScalarParameterLibrary>& paramLib) : 
      alpha(factor) 
    {
      // Add nonlinear factor to parameter library
      std::string name = "Exponential Source Function Nonlinear Factor";
      if (!paramLib->isParameter(name))
	paramLib->addParameterFamily(name, true, false);
      if (!paramLib->template isParameterForType<ScalarT>(name)) {
	Teuchos::RCP< ExponentialNonlinearFactorParameter<ScalarT> > tmp = 
	  Teuchos::rcp(new ExponentialNonlinearFactorParameter<ScalarT>(Teuchos::rcp(this,false)));
	paramLib->template addEntry<ScalarT>(name, tmp);
      }
    };

    //! Destructor
    virtual ~ExponentialSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
	     std::vector<ScalarT>& value) const {
      for (unsigned int i=0; i<solution.size(); i++) {
	value[i] = -std::exp(alpha)*std::exp(solution[i]);
	//value[i] = -1.0;
      }
      
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
    ExponentialSourceFunction(const ExponentialSourceFunction&);

    //! Private to prohibit copying
    ExponentialSourceFunction& operator=(const ExponentialSourceFunction&);

  protected:
  
    //! Factor
    ScalarT alpha;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * the nonlinear factor in the cubic source function
   */
  template <typename ScalarT>
  class ExponentialNonlinearFactorParameter : 
    public Sacado::ScalarParameterEntry<ScalarT> {

  public:

    //! Constructor
    ExponentialNonlinearFactorParameter(
			const Teuchos::RCP< ExponentialSourceFunction<ScalarT> >& s) 
      : srcFunc(s) {}

    //! Destructor
    virtual ~ExponentialNonlinearFactorParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      setValueAsConstant(ScalarT(value)); }

    //! Set parameter this object represents to \em value
    virtual void setValueAsConstant(const ScalarT& value) { 
      srcFunc->setFactor(value, true); }
    
    //! Set parameter this object represents to \em value
    virtual void setValueAsIndependent(const ScalarT& value) { 
      srcFunc->setFactor(value, false); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { return srcFunc->getFactor(); }
    
  protected:  
    
    //! Pointer to source function
    Teuchos::RCP< ExponentialSourceFunction<ScalarT> > srcFunc;

  };

}

#endif // FEAPP_CUBICSOURCEFUNCTION_HPP
