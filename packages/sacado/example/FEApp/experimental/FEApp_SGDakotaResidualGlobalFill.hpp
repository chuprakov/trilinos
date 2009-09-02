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

#ifndef FEAPP_SGDAKOTARESIDUALGLOBALFILL_HPP
#define FEAPP_SGDAKOTARESIDUALGLOBALFILL_HPP

#include "FEApp_TemplateTypes.hpp"
#if SG_ACTIVE

#include "Stokhos_ConfigDefs.h"
#ifdef HAVE_DAKOTA

#include "FEApp_GlobalFill.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Sacado_ScalarParameterVector.hpp"

// Dakota includes
#include "system_defs.h"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "FEApp_DakotaElementResidualInterface.hpp"

namespace FEApp {

  class SGDakotaResidualGlobalFill : public GlobalFill<SGResidualType> {
  public:

    //! Scalar type
    typedef FEApp::EvaluationTraits::apply<SGResidualType>::type ScalarT;
    
    //! Constructor
    SGDakotaResidualGlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<SGResidualType> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sgBasis,
      const Teuchos::RCP< FEApp::AbstractPDE<ResidualType> >& resPDEEquations,
      const ParamVec* pvec);
  
    //! Destructor
    virtual ~SGDakotaResidualGlobalFill();

    //! Compute global fill
    virtual void 
    computeGlobalFill(FEApp::AbstractInitPostOp<SGResidualType>& initPostOp);

  private:

    //! Private to prohibit copying
    SGDakotaResidualGlobalFill(const SGDakotaResidualGlobalFill&);

    //! Private to prohibit copying
    SGDakotaResidualGlobalFill& operator=(const SGDakotaResidualGlobalFill&);

  protected:
    
    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;
    Teuchos::RCP< FEApp::AbstractPDE<ResidualType> > residPDE;
    const ParamVec* p;
    unsigned int sg_size;
    Dakota::ParallelLibrary parallel_lib;
    Dakota::ProblemDescDB problem_db;
    Dakota::Strategy selected_strategy;
    FEApp::DakotaElementResidualInterface *elemInterface;
  };

}

#endif // HAVE_DAKOTA

#endif // SG_ACTIVE

#endif // SGDAKOTARESIDUALGLOBALFILL_HPP
