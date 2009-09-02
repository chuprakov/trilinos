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

#ifndef FEAPP_BRUSSELATORPDE_HPP
#define FEAPP_BRUSSELATORPDE_HPP

#include "Teuchos_RCP.hpp"

#include "FEApp_AbstractPDE.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"

namespace FEApp {

  template <typename EvalT>
  class BrusselatorPDE : public FEApp::AbstractPDE<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractPDE<EvalT>::ScalarT ScalarT;
  
    //! Constructor
    BrusselatorPDE(
                double alpha, double beta, double D1, double D2,
                const Teuchos::RCP<ParamLib>& paramLib);

    //! Destructor
    virtual ~BrusselatorPDE();

    //! Number of discretized equations
    virtual unsigned int numEquations() const;

    //! Initialize PDE
    virtual void init(unsigned int numQuadPoints, unsigned int numNodes);

    //! Evaluate discretized PDE element-level residual
    virtual void
    evaluateElementResidual(const FEApp::AbstractQuadrature& quadRule,
                            const FEApp::AbstractElement& element,
                            const std::vector<ScalarT>* dot,
                            const std::vector<ScalarT>& solution,
                            std::vector<ScalarT>& residual);

  private:
    
    //! Private to prohibit copying
    BrusselatorPDE(const BrusselatorPDE&);
    
    //! Private to prohibit copying
    BrusselatorPDE& operator=(const BrusselatorPDE&);
    
  protected:
    
    //! Number of quad points
    unsigned int num_qp;
    
    //! Number of nodes
    unsigned int num_nodes;
    
    //! Shape function values
    std::vector< std::vector<double> > phi;

    //! Shape function derivatives
    std::vector< std::vector<double> > dphi;

    //! Element transformation Jacobian
    std::vector<double> jac;

    //! Discretized solution
    std::vector<ScalarT> T;

    //! Discretized solution
    std::vector<ScalarT> C;

    //! Discretized solution derivative
    std::vector<ScalarT> dT;

    //! Discretized solution derivative
    std::vector<ScalarT> dC;

    //! Discretized time derivative
    std::vector<ScalarT> Tdot;

    //! Discretized time derivative
    std::vector<ScalarT> Cdot;

    //! Model parameters
    ScalarT alpha, beta;
    double D1, D2;

    //! Parameter library
    Teuchos::RCP<ParamLib> pl;

  };

  class BrusselatorPDE_TemplateBuilder {
  public:
    BrusselatorPDE_TemplateBuilder(
          double alpha_, double beta_, double D1_, double D2_,
          const Teuchos::RCP<ParamLib>& paramLib) :
      alpha(alpha_), beta(beta_), D1(D1_), D2(D2_), pl(paramLib) {}
    template <typename T>
    Teuchos::RCP<FEApp::AbstractPDE_NTBase> build() const {
      return Teuchos::rcp( new BrusselatorPDE<T>(alpha, beta, D1, D2, pl));
    }
  protected:
    double alpha, beta, D1, D2;
    Teuchos::RCP<ParamLib> pl;
  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_BrusselatorPDEImpl.hpp"
#endif 

#endif // FEAPP_BRUSSELATORPDE_HPP
