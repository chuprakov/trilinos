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

#ifndef FEAPP_BRUSSELATORPROBLEM_HPP
#define FEAPP_BRUSSELATORPROBLEM_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "FEApp_AbstractProblem.hpp"
#include "FEApp_BrusselatorPDE.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"

namespace FEApp {

  /*!
   * \brief Abstract interface for representing a 1-D finite element
   * problem.
   */
  class BrusselatorProblem : public FEApp::AbstractProblem {
  public:
  
    //! Default constructor
    BrusselatorProblem(
              const Teuchos::RCP<Teuchos::ParameterList>& params,
              const Teuchos::RCP<ParamLib>& paramLib);

    //! Destructor
    virtual ~BrusselatorProblem();

    //! Get the number of equations
    virtual unsigned int numEquations() const;

    //! Build the PDE instantiations, boundary conditions, and initial solution
    virtual void 
    buildProblem(const Epetra_Map& dofMap,
                 const Epetra_Map& overlapped_dofMap,
                 FEApp::AbstractPDE_TemplateManager<EvalTypes>& pdeTM,
                 std::vector< Teuchos::RCP<FEApp::NodeBC> >& bcs,
                 const Teuchos::RCP<Epetra_Vector>& u);

  private:

    //! Private to prohibit copying
    BrusselatorProblem(const BrusselatorProblem&);
    
    //! Private to prohibit copying
    BrusselatorProblem& operator=(const BrusselatorProblem&);

  protected:

    //! Model parameters
    double alpha, beta, D1, D2;

    //! Parameter library
    Teuchos::RCP<ParamLib> paramLib;

  };

  

}

#endif // FEAPP_BRUSSELATORPROBLEM_HPP
