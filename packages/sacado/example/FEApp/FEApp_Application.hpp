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

#ifndef FEAPP_APPLICATION_HPP
#define FEAPP_APPLICATION_HPP

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"

#include "FEApp_NodeBC.hpp"
#include "FEApp_AbstractPDE.hpp"
#include "FEApp_AbstractQuadrature.hpp"
#include "FEApp_AbstractDiscretization.hpp"
#include "FEApp_AbstractProblem.hpp"
#include "FEApp_TemplateTypes.hpp"
#include "FEApp_InitPostOps.hpp"

#include "Sacado_ScalarParameterLibrary.hpp"
#include "Sacado_ScalarParameterVector.hpp"

#if SG_ACTIVE
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "FEApp_BlockDiscretization.hpp"
#endif

#if SGFAD_ACTIVE
#include "Stokhos_TripleProduct.hpp"
#endif

namespace FEApp {

  class Application {
  public:

    //! Constructor 
    Application(const std::vector<double>& coords,
		const Teuchos::RCP<const Epetra_Comm>& comm,
		const Teuchos::RCP<Teuchos::ParameterList>& params,
		bool is_transient);

    //! Destructor
    ~Application();

    //! Get DOF map
    Teuchos::RCP<const Epetra_Map> getMap() const;

    //! Get Jacobian graph
    Teuchos::RCP<const Epetra_CrsGraph> getJacobianGraph() const;

    //! Get initial solution
    Teuchos::RCP<const Epetra_Vector> getInitialSolution() const;

    //! Get parameter library
    Teuchos::RCP<Sacado::ScalarParameterLibrary> getParamLib();

    //! Return whether problem is transient
    bool isTransient() const;

    //! Create new W operator
    Teuchos::RCP<Epetra_Operator> createW() const;

    //! Create new preconditioner operator
    Teuchos::RCP<Epetra_Operator> createPrec() const;

    //! Compute global residual
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalResidual(const Epetra_Vector* xdot,
			       const Epetra_Vector& x,
			       const Sacado::ScalarParameterVector* p,
			       Epetra_Vector& f);

    //! Compute global Jacobian
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalJacobian(double alpha, double beta,
			       const Epetra_Vector* xdot,
			       const Epetra_Vector& x,
			       const Sacado::ScalarParameterVector* p,
			       Epetra_Vector* f,
			       Epetra_Operator& jac);

    //! Compute global Preconditioner
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalPreconditioner(double alpha, double beta,
				     const Epetra_Vector* xdot,
				     const Epetra_Vector& x,
				     const Sacado::ScalarParameterVector* p,
				     Epetra_Vector* f,
				     Epetra_Operator& jac);
    
    //! Compute global Tangent
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalTangent(double alpha, double beta,
			      bool sum_derivs,
			      const Epetra_Vector* xdot,
			      const Epetra_Vector& x,
			      Sacado::ScalarParameterVector* p,
			      const Epetra_MultiVector* Vx,
			      const Teuchos::SerialDenseMatrix<int,double>* Vp,
			      Epetra_Vector* f,
			      Epetra_MultiVector* JVx,
			      Epetra_MultiVector* fVp);

    //! Compute global residual for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGResidual(const Epetra_Vector* sg_xdot,
				 const Epetra_Vector& sg_x,
				 const Sacado::ScalarParameterVector* p,
				 Epetra_Vector& sg_f);

    //! Compute global Jacobian for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGJacobian(double alpha, double beta,
				 const Epetra_Vector* sg_xdot,
				 const Epetra_Vector& sg_x,
				 const Sacado::ScalarParameterVector* p,
				 Epetra_Vector* sg_f,
				 Epetra_Operator& sg_jac); 

    //! Compute global Preconditioner for stochastic Galerkin problem
    /*!
     * Set xdot to NULL for steady-state problems
     */
    void computeGlobalSGPreconditioner(double alpha, double beta,
				       const Epetra_Vector* sg_xdot,
				       const Epetra_Vector& sg_x,
				       const Sacado::ScalarParameterVector* p,
				       Epetra_Vector* sg_f,
				       Epetra_Operator& sg_jac); 

  private:
    
    //! Private to prohibit copying
    Application(const Application&);

    //! Private to prohibit copying
    Application& operator=(const Application&);

  protected:

    //! Is problem transient
    bool transient;
    
    //! Element discretization
    Teuchos::RCP<FEApp::AbstractDiscretization> disc;

    //! Boundary conditions
    std::vector< Teuchos::RCP<FEApp::NodeBC> > bc;

    //! Quadrature rule
    Teuchos::RCP<const FEApp::AbstractQuadrature> quad;

    //! PDE equations
    FEApp::AbstractPDE_TemplateManager<ValidTypes> pdeTM;

    //! Initial solution vector
    Teuchos::RCP<Epetra_Vector> initial_x;

    //! Importer for overlapped data
    Teuchos::RCP<Epetra_Import> importer;

    //! Exporter for overlapped data
    Teuchos::RCP<Epetra_Export> exporter;

    //! Overlapped solution vector
    Teuchos::RCP<Epetra_Vector> overlapped_x;

    //! Overlapped time derivative vector
    Teuchos::RCP<Epetra_Vector> overlapped_xdot;

    //! Overlapped residual vector
    Teuchos::RCP<Epetra_Vector> overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> overlapped_jac;

    //! Parameter library
    Teuchos::RCP<Sacado::ScalarParameterLibrary> paramLib;

    //! Enable stochastic Galerkin discretization
    bool enable_sg;

#if SG_ACTIVE

    //! Solver method
    std::string sg_solver_method;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > sg_basis;

     typedef SGType::expansion_type::tp_type tp_type;
    //! Stochastic Galerkin triple product
    Teuchos::RCP<const tp_type> Cijk;

    //! Stochastic Galerkin discretization
    Teuchos::RCP<FEApp::BlockDiscretization> sg_disc;

    //! Initial solution vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_initial_x;

    //! Importer for overlapped data
    Teuchos::RCP<Epetra_Import> sg_importer;

    //! Exporter for overlapped data
    Teuchos::RCP<Epetra_Export> sg_exporter;

    //! Overlapped solution vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_overlapped_x;

    //! Overlapped time derivative vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_overlapped_xdot;

    //! Overlapped residual vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_overlapped_f;

    //! Overlapped Jacobian matrix
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> sg_overlapped_jac;

    //! SG Residual fill op
    Teuchos::RCP<FEApp::SGResidualOp> sg_res_fill_op;

    //! SG Jacobian fill op
    Teuchos::RCP<FEApp::SGJacobianOp> sg_full_jac_fill_op;

    //! SG Matrix Free Jacobian fill op
    Teuchos::RCP<FEApp::SGMatrixFreeJacobianOp> sg_mf_jac_fill_op;

    //! SG Preconditioner parameters
    Teuchos::RCP<Teuchos::ParameterList> precParams;

#endif

  };

}

#endif // FEAPP_APPLICATION_HPP
