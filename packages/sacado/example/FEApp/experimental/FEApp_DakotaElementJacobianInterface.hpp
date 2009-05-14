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

#ifndef FEAPP_DAKOTAELEMENTJACOBIANINTERFACE_HPP
#define FEAPP_DAKOTAELEMENTJACOBIANINTERFACE_HPP

#include "Sacado_ConfigDefs.h"
#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_DAKOTA

#include "DirectApplicInterface.H"
#include "CommandLineHandler.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"

#include <vector>
#include "Teuchos_RCP.hpp"
#include "FEApp_TemplateTypes.hpp"
#include "FEApp_AbstractPDE.hpp"
#include "FEApp_AbstractQuadrature.hpp"
#include "FEApp_AbstractElement.hpp"

namespace FEApp {

  class DakotaElementJacobianInterface : public Dakota::DirectApplicInterface {
  public:

    //! Constructor
    DakotaElementJacobianInterface(
         const Dakota::ProblemDescDB& problem_db_,
			   const Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> >& pde_,
         const ParamVec& pvec_,
         const Teuchos::RCP<const FEApp::AbstractQuadrature>& quad_,
         unsigned int ndof_,
         double alpha_,
         double beta_);

    //! Destructor
    ~DakotaElementJacobianInterface() {};

    //! Reset interface for a new element
    void reset(const Teuchos::RCP<const FEApp::AbstractElement>& e_,
               const Teuchos::RCP< std::vector<SGFadType> >& sg_xdot_,
               const Teuchos::RCP< std::vector<SGFadType> >& sg_x_);

  protected:

    // Virtual function redefinitions

    //int derived_map_if(const Dakota::String& if_name);
    int derived_map_ac(const Dakota::String& ac_name);
    //int derived_map_of(const Dakota::String& of_name);

  private:

    // Data
    Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > pde;
    ParamVec pvec;
    Teuchos::RCP<const FEApp::AbstractQuadrature> quad;
    unsigned int numParameters;
    unsigned int numEquations;
    Teuchos::RCP<const FEApp::AbstractElement> e;
    Teuchos::RCP< std::vector<SGFadType> > sg_xdot;
    Teuchos::RCP< std::vector<SGFadType> > sg_x;
    std::vector<FadType> xdot;
    std::vector<FadType> x;
    std::vector<FadType> f;
    double alpha;
    double beta;
  };

} // namespace FEApp

#endif // HAVE_DAKOTA

#endif // FEAPP_DAKOTAELEMENTJACOBIANLINTERFACE_HPP
