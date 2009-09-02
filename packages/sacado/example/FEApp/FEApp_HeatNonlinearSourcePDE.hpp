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

#ifndef FEAPP_HEATNONLINEARSOURCEPDE_HPP
#define FEAPP_HEATNONLINEARSOURCEPDE_HPP

#include "Teuchos_RCP.hpp"

#include "FEApp_AbstractPDE.hpp"
#include "FEApp_FunctionFactory.hpp"
#include "FEApp_SourceFunctionFactory.hpp"
#include "Sacado_ScalarParameterLibrary.hpp"

namespace FEApp {

  template <typename EvalT>
  class HeatNonlinearSourcePDE : public FEApp::AbstractPDE<EvalT> {
  public:

    //! Scalar type
    typedef typename FEApp::AbstractPDE<EvalT>::ScalarT ScalarT;
  
    //! Constructor
    HeatNonlinearSourcePDE(
      const Teuchos::RCP< const FEApp::AbstractFunction<EvalT> >& mat_func,
			const Teuchos::RCP< const FEApp::AbstractSourceFunction<EvalT> >& src_func);

    //! Destructor
    virtual ~HeatNonlinearSourcePDE();

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
    HeatNonlinearSourcePDE(const HeatNonlinearSourcePDE&);

    //! Private to prohibit copying
    HeatNonlinearSourcePDE& operator=(const HeatNonlinearSourcePDE&);

  protected:

    //! Pointer to material function
    Teuchos::RCP< const FEApp::AbstractFunction<EvalT> > mat;
    
    //! Pointer to source function
    Teuchos::RCP< const FEApp::AbstractSourceFunction<EvalT> > source;

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
    std::vector<ScalarT> u;

    //! Discretized solution
    std::vector<ScalarT> du;

    //! Discretized time derivative
    std::vector<ScalarT> udot;

    //! Material function values
    std::vector<ScalarT> a;

    //! Source function values
    std::vector<ScalarT> f;

  };

  class HeatNonlinearSourcePDE_TemplateBuilder {
  public:
    HeatNonlinearSourcePDE_TemplateBuilder(
		const Teuchos::RCP<Teuchos::ParameterList>& params_,
	        const Teuchos::RCP<ParamLib>& paramLib) :
      mat_params(Teuchos::rcp(&(params_->sublist("Material Function")),false)),
      src_params(Teuchos::rcp(&(params_->sublist("Source Function")),false)),
      pl(paramLib) {}
    template <typename T>
    Teuchos::RCP<FEApp::AbstractPDE_NTBase> build() const {
      FEApp::FunctionFactory<T> matFactory(mat_params, pl);
      Teuchos::RCP< FEApp::AbstractFunction<T> > mat =
        matFactory.create();
      FEApp::SourceFunctionFactory<T> srcFactory(src_params, pl);
      Teuchos::RCP< FEApp::AbstractSourceFunction<T> > source =
        srcFactory.create();
      return Teuchos::rcp( new FEApp::HeatNonlinearSourcePDE<T>(mat, source));
    }
  protected:
    Teuchos::RCP<Teuchos::ParameterList> mat_params;
    Teuchos::RCP<Teuchos::ParameterList> src_params;
    Teuchos::RCP<ParamLib> pl;
  };

}

// Include implementation
#ifndef SACADO_ETI
#include "FEApp_HeatNonlinearSourcePDEImpl.hpp"
#endif 

#endif // FEAPP_HEATNONLINERASOURCEPDE_HPP
