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

#ifndef FEAPP_MODELEVALUATOR_HPP
#define FEAPP_MODELEVALUATOR_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"

#include "EpetraExt_ModelEvaluator.h"

#include "FEApp_Application.hpp"

#include "Epetra_LocalMap.h"
#include "Sacado_ScalarParameterVector.hpp"

namespace FEApp {

  class ModelEvaluator : public EpetraExt::ModelEvaluator {
  public:

    // Constructor
    ModelEvaluator(
       const Teuchos::RCP<FEApp::Application>& app,
       const Teuchos::RCP< Teuchos::Array<std::string> >& free_param_names = 
       Teuchos::null,
       const Teuchos::RCP< Teuchos::Array<std::string> >& sg_param_names = 
       Teuchos::null);

    /** \name Overridden from EpetraExt::ModelEvaluator . */
    //@{

    //! Return solution vector map
    Teuchos::RCP<const Epetra_Map> get_x_map() const;

    //! Return residual vector map
    Teuchos::RCP<const Epetra_Map> get_f_map() const;

    //! Return parameter vector map
    Teuchos::RCP<const Epetra_Map> get_p_map(int l) const;

    //! Return array of parameter names
    Teuchos::RCP<const Teuchos::Array<std::string> > 
    get_p_names(int l) const;

    //! Return initial solution
    Teuchos::RCP<const Epetra_Vector> get_x_init() const;

    //! Return initial parameters
    Teuchos::RCP<const Epetra_Vector> get_p_init(int l) const;

    //! Create W = alpha*M + beta*J matrix
    Teuchos::RCP<Epetra_Operator> create_W() const;

    //! Create InArgs
    InArgs createInArgs() const;

    //! Create OutArgs
    OutArgs createOutArgs() const;

    //! Evaluate model on InArgs
    void evalModel(const InArgs& inArgs, const OutArgs& outArgs) const;

    //@}

  protected:

    //! Application object
    Teuchos::RCP<FEApp::Application> app;

    //! List of free parameter names
    Teuchos::Array< Teuchos::RCP< Teuchos::Array<std::string> > > param_names;

    //! Sacado parameter vector
    Teuchos::Array< Teuchos::RCP<ParamVec> > sacado_param_vec;

    //! Epetra map for parameter vector
    Teuchos::Array< Teuchos::RCP<Epetra_LocalMap> > epetra_param_map;

    //! Epetra parameter vector
    Teuchos::Array< Teuchos::RCP<Epetra_Vector> > epetra_param_vec;

    //! Whether we are support SG
    bool supports_sg;

#if SG_ACTIVE
    //! Stochastic Galerkin parameters
    mutable Teuchos::Array<SGType> p_sg_vals;
#endif

  };

}

#endif // FEAPP_MODELEVALUATOR_HPP
