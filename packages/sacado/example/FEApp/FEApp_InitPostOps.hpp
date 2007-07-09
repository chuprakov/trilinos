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

#ifndef FEAPP_INITPOSTOPS_HPP
#define FEAPP_INITPOSTOPS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "FEApp_AbstractInitPostOp.hpp"

#include "Sacado_Fad_DFad.hpp"
#include "Sacado_ScalarParameterVector.hpp"

namespace FEApp {

  //! Fill operator for residual
  class ResidualOp : public FEApp::AbstractInitPostOp<double> {
  public:
    
    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    ResidualOp(
	    const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
	    const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
	    const Teuchos::RCP<Epetra_Vector>& overlapped_f);

    //! Destructor
    virtual ~ResidualOp();
    
    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
			     unsigned int neqn,
			     std::vector<double>* elem_xdot,
			     std::vector<double>& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
			     unsigned int neqn,
			     std::vector<double>& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector<double>* node_xdot,
			  std::vector<double>& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector<double>& node_f);

  private:
    
    //! Private to prohibit copying
    ResidualOp(const ResidualOp&);

    //! Private to prohibit copying
    ResidualOp& operator=(const ResidualOp&);

  protected:

    //! Time derivative vector (may be null)
    Teuchos::RCP<const Epetra_Vector> xdot;

    //! Solution vector
    Teuchos::RCP<const Epetra_Vector> x;

    //! Residual vector
    Teuchos::RCP<Epetra_Vector> f;

  };

  //! Fill operator for Jacobian
  class JacobianOp : 
    public FEApp::AbstractInitPostOp< Sacado::Fad::DFad<double> > {
  public:

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    JacobianOp(
	    double alpha, double beta,
	    const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
	    const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
	    const Teuchos::RCP<Epetra_Vector>& overlapped_f,
	    const Teuchos::RCP<Epetra_CrsMatrix>& overlapped_jac);

    //! Destructor
    virtual ~JacobianOp();

    //! Evaulate element init operator
    virtual void elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< Sacado::Fad::DFad<double> >* elem_xdot,
			 std::vector< Sacado::Fad::DFad<double> >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
			     unsigned int neqn,
			     std::vector< Sacado::Fad::DFad<double> >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector< Sacado::Fad::DFad<double> >* node_xdot,
			  std::vector< Sacado::Fad::DFad<double> >& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector< Sacado::Fad::DFad<double> >& node_f);

  private:
    
    //! Private to prohibit copying
    JacobianOp(const JacobianOp&);
    
    //! Private to prohibit copying
    JacobianOp& operator=(const JacobianOp&);

  protected:

    //! Coefficient of mass matrix
    double m_coeff;

    //! Coefficient of Jacobian matrix
    double j_coeff;

    //! Time derivative vector (may be null)
    Teuchos::RCP<const Epetra_Vector> xdot;

    //! Solution vector
    Teuchos::RCP<const Epetra_Vector> x;

    //! Residual vector
    Teuchos::RCP<Epetra_Vector> f;

    //! Jacobian matrix
    Teuchos::RCP<Epetra_CrsMatrix> jac;

  };

  /*!! 
   * Fill operator for "tangent":  
   * alpha*df/dxdot*Vxdot + beta*df/dx*Vx + df/dp*V_p
   */
  class TangentOp : 
    public FEApp::AbstractInitPostOp< Sacado::Fad::DFad<double> > {
  public:

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    TangentOp(
	double alpha, double beta, bool sum_derivs,
	const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
	const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
	const Teuchos::RCP<Sacado::ScalarParameterVector>& p,
	const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vx,
	const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vxdot,
	const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp,
	const Teuchos::RCP<Epetra_Vector>& overlapped_f,
	const Teuchos::RCP<Epetra_MultiVector>& overlapped_JV,
	const Teuchos::RCP<Epetra_MultiVector>& overlapped_fp);

    //! Destructor
    virtual ~TangentOp();

    //! Evaulate element init operator
    virtual void elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< Sacado::Fad::DFad<double> >* elem_xdot,
			 std::vector< Sacado::Fad::DFad<double> >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
			     unsigned int neqn,
			     std::vector< Sacado::Fad::DFad<double> >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector< Sacado::Fad::DFad<double> >* node_xdot,
			  std::vector< Sacado::Fad::DFad<double> >& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
			  unsigned int neqn,
			  std::vector< Sacado::Fad::DFad<double> >& node_f);

  private:
    
    //! Private to prohibit copying
    TangentOp(const TangentOp&);
    
    //! Private to prohibit copying
    TangentOp& operator=(const TangentOp&);

  protected:

    //! Coefficient of mass matrix
    double m_coeff;

    //! Coefficient of Jacobian matrix
    double j_coeff;

    //! Whether to sum derivative terms
    bool sum_derivs;

    //! Time derivative vector (may be null)
    Teuchos::RCP<const Epetra_Vector> xdot;

    //! Solution vector
    Teuchos::RCP<const Epetra_Vector> x;

    //! Parameter vector for parameter derivatives
    Teuchos::RCP<Sacado::ScalarParameterVector> params;

    //! Seed matrix for state variables
    Teuchos::RCP<const Epetra_MultiVector> Vx;

    //! Seed matrix for transient variables
    Teuchos::RCP<const Epetra_MultiVector> Vxdot;

    //! Seed matrix for parameters
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> > Vp;

    //! Residual vector
    Teuchos::RCP<Epetra_Vector> f;

    //! Tangent matrix (alpha*df/dxdot + beta*df/dx)*V
    Teuchos::RCP<Epetra_MultiVector> JV;
    
    //! Tangent matrix df/dp*V_p
    Teuchos::RCP<Epetra_MultiVector> fp;

    //! Stores number of columns in seed matrix V
    int num_cols_x;

    //! Stores number of columns in seend matrix Vp
    int num_cols_p;

    //! Stores the total number of columns
    int num_cols_tot;

    //! Stores the parameter offset
    int param_offset;

  };

}

#endif // FEAPP_INITPOSTOPS_HPP
