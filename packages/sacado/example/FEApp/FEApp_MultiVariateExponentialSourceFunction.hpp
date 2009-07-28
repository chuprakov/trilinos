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

#ifndef FEAPP_MULTIVARIATEEXPONENTIALSOURCEFUNCTION_HPP
#define FEAPP_MULTIVARIATEEXPONENTIALSOURCEFUNCTION_HPP

#include <sstream>

#include "FEApp_AbstractSourceFunction.hpp"

#include "Teuchos_RCP.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterEntry.hpp"
#include "Sacado_Traits.hpp"

namespace FEApp {

  template <typename EvalT> 
  class MultiVariateExponentialNonlinearFactorParameter;

  /*!
   * \brief An exponential PDE source function
   */
  template <typename EvalT>
  class MultiVariateExponentialSourceFunction : 
    public FEApp::AbstractSourceFunction<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractSourceFunction<EvalT>::ScalarT ScalarT;
  
    //! Default constructor
    MultiVariateExponentialSourceFunction(
         const std::vector<ScalarT>& factor,
	       const Teuchos::RCP<ParamLib>& paramLib) : 
      alpha(factor) 
    {
      // Add nonlinear factor to parameter library
      for (unsigned int i=0; i<factor.size(); i++) {
        std::stringstream ss;
        ss << "Exponential Source Function Nonlinear Factor " << i;
        std::string name = ss.str();
        if (!paramLib->isParameter(name))
          paramLib->addParameterFamily(name, true, false);
        if (!paramLib->template isParameterForType<EvalT>(name)) {
          Teuchos::RCP< MultiVariateExponentialNonlinearFactorParameter<EvalT> > tmp = 
            Teuchos::rcp(new MultiVariateExponentialNonlinearFactorParameter<EvalT>(Teuchos::rcp(this,false),i));
          paramLib->template addEntry<EvalT>(name, tmp);
        }
      }
    };

    //! Destructor
    virtual ~MultiVariateExponentialSourceFunction() {};

    //! Evaluate source function
    virtual void
    evaluate(const std::vector<ScalarT>& solution,
             std::vector<ScalarT>& value) const {
      ScalarT a = 0.0;
      for (unsigned int j=0; j<alpha.size(); j++)
	a += alpha[j];
      a /= static_cast<double>(alpha.size());
      for (unsigned int i=0; i<solution.size(); i++) {
        value[i] = a*std::exp(solution[i]);
      }
      
    }

    //! Set nonlinear factor
    void setFactor(unsigned int i, const ScalarT& val, bool mark_constant) { 
      alpha[i] = val;
      if (mark_constant) Sacado::MarkConstant<ScalarT>::eval(alpha[i]); 
    }

    //! Get nonlinear factor
    const ScalarT& getFactor(unsigned int i) const { return alpha[i]; }

  private:

    //! Private to prohibit copying
    MultiVariateExponentialSourceFunction(const MultiVariateExponentialSourceFunction&);

    //! Private to prohibit copying
    MultiVariateExponentialSourceFunction& operator=(const MultiVariateExponentialSourceFunction&);

  protected:
  
    //! Factor
    std::vector<ScalarT> alpha;

  };

  /*!
   * @brief Parameter class for sensitivity/stability analysis representing
   * the nonlinear factor in the cubic source function
   */
  template <typename EvalT>
  class MultiVariateExponentialNonlinearFactorParameter : 
    public Sacado::ScalarParameterEntry<EvalT,EvaluationTraits> {

  public:

    //! Scalar type
    typedef typename Sacado::ScalarParameterEntry<EvalT,EvaluationTraits>::ScalarT ScalarT;

    //! Constructor
    MultiVariateExponentialNonlinearFactorParameter(
        const Teuchos::RCP< MultiVariateExponentialSourceFunction<EvalT> >& s,
        unsigned int i) 
      : srcFunc(s), index(i) {}

    //! Destructor
    virtual ~MultiVariateExponentialNonlinearFactorParameter() {}

    //! Set real parameter value
    virtual void setRealValue(double value) { 
      srcFunc->setFactor(index, value, true); }
    
    //! Set parameter this object represents to \em value
    virtual void setValue(const ScalarT& value) { 
      srcFunc->setFactor(index, value, false); }

    //! Get real parameter value
    virtual double getRealValue() const {
      return Sacado::ScalarValue<ScalarT>::eval(srcFunc->getFactor(index)); }
    
    //! Get parameter value this object represents
    virtual const ScalarT& getValue() const { 
      return srcFunc->getFactor(index); }
    
  protected:  
    
    //! Pointer to source function
    Teuchos::RCP< MultiVariateExponentialSourceFunction<EvalT> > srcFunc;

    //! Index of factor
    unsigned int index;

  };

}

#endif // FEAPP_MULTIVARIATEEXPONENTIALSOURCEFUNCTION_HPP
