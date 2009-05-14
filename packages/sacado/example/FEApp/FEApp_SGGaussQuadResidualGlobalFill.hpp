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

#ifndef FEAPP_SGGAUSSQUADRESIDUALGLOBALFILL_HPP
#define FEAPP_SGGAUSSQUADRESIDUALGLOBALFILL_HPP

#include "FEApp_TemplateTypes.hpp"
#if SG_ACTIVE

#include "FEApp_GlobalFill.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"
#include "Sacado_ScalarParameterVector.hpp"
#include "Teuchos_BLAS.hpp"

namespace FEApp {

  class SGGaussQuadResidualGlobalFill : public GlobalFill<SGResidualType> {
  public:

    //! Scalar type
    typedef FEApp::EvaluationTraits::apply<SGResidualType>::type ScalarT;
    
    //! Constructor
    SGGaussQuadResidualGlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<SGResidualType> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sgBasis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sgQuad,
      const Teuchos::RCP< FEApp::AbstractPDE<ResidualType> >& resPDEEquations,
      const ParamVec* pvec);
  
    //! Destructor
    virtual ~SGGaussQuadResidualGlobalFill();

    //! Compute global fill
    virtual void 
    computeGlobalFill(FEApp::AbstractInitPostOp<SGResidualType>& initPostOp);

  private:

    //! Private to prohibit copying
    SGGaussQuadResidualGlobalFill(const SGGaussQuadResidualGlobalFill&);

    //! Private to prohibit copying
    SGGaussQuadResidualGlobalFill& operator=(const SGGaussQuadResidualGlobalFill&);

  protected:
    
    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Stochastic Galerkin quadrature
    Teuchos::RCP<const Stokhos::Quadrature<int,double> > sg_quad;

    Teuchos::RCP< FEApp::AbstractPDE<ResidualType> > residPDE;
    const ParamVec* p;
    const std::vector< std::vector<double> >& quad_points;
    const std::vector<double>& quad_weights;
    const std::vector< std::vector<double> >& quad_values;
    const std::vector<double>& norms;
    unsigned int sg_size;
    unsigned int nqp;
    std::vector<double> x;
    std::vector<double>* xdot;
    std::vector<double> f;

    std::vector<double> xqp;
    std::vector<double> xdotqp;
    std::vector<double> pqp;
    std::vector<double> fqp;

    std::vector<double> qv;
    std::vector<double> sqv;

    std::vector<double> sg_x;
    std::vector<double> sg_xdot;
    std::vector<double> sg_p;
    std::vector<double> sg_f;

    //! BLAS wrappers
    Teuchos::BLAS<int,double> blas;

  };

}

#endif // SG_ACTIVE

#endif // SGGAUSSQUADRESIDUALGLOBALFILL_HPP
