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
#include "FEApp_TemplateTypes.hpp"
#include "Sacado_ScalarParameterVector.hpp"

#if SG_ACTIVE
#include "EpetraExt_BlockVector.h"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "EpetraExt_BlockCrsMatrix.h"
#endif

namespace FEApp {

  //! Fill operator for residual
  class ResidualOp : public FEApp::AbstractInitPostOp<FEApp::ResidualType> {
  public:
    
    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    ResidualOp(const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
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

    //! Finalize fill
    virtual void finalizeFill() {}

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
    public FEApp::AbstractInitPostOp<FEApp::JacobianType> {
  public:

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    JacobianOp(double alpha, double beta,
               const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
               const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
               const Teuchos::RCP<Epetra_Vector>& overlapped_f,
               const Teuchos::RCP<Epetra_CrsMatrix>& overlapped_jac);

    //! Destructor
    virtual ~JacobianOp();

    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< FadType >* elem_xdot,
                             std::vector< FadType >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< FadType >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< FadType >* node_xdot,
                          std::vector< FadType >& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< FadType >& node_f);

    //! Finalize fill
    virtual void finalizeFill() {}

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
    public FEApp::AbstractInitPostOp<FEApp::TangentType> {
  public:

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    TangentOp(
	       double alpha, double beta, bool sum_derivs,
         const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
         const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
         const Teuchos::RCP<ParamVec>& p,
         const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vx,
         const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vxdot,
	       const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp,
         const Teuchos::RCP<Epetra_Vector>& overlapped_f,
         const Teuchos::RCP<Epetra_MultiVector>& overlapped_JV,
         const Teuchos::RCP<Epetra_MultiVector>& overlapped_fp);

    //! Destructor
    virtual ~TangentOp();

    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< FadType >* elem_xdot,
                             std::vector< FadType >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< FadType >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< FadType >* node_xdot,
                          std::vector< FadType >& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< FadType >& node_f);

    //! Finalize fill
    virtual void finalizeFill() {}

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
    Teuchos::RCP<ParamVec> params;

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

#if SG_ACTIVE

  //! Fill operator for Stochastic Galerkin residual
  class SGResidualOp : 
    public FEApp::AbstractInitPostOp<FEApp::SGResidualType> {
  public:
    
    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    SGResidualOp(
	      const Teuchos::RCP<const Epetra_Map>& base_map,
	      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis,
        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
        const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f);

    //! Destructor
    virtual ~SGResidualOp();

    //! Reset operator for new fill
    virtual void reset(
	        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
          const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x);
    
    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector<SGType>* elem_xdot,
                             std::vector<SGType>& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector<SGType>& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector<SGType>* node_xdot,
                          std::vector<SGType>& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector<SGType>& node_f);

    //! Finalize fill
    virtual void finalizeFill();

  private:
    
    //! Private to prohibit copying
    SGResidualOp(const SGResidualOp&);

    //! Private to prohibit copying
    SGResidualOp& operator=(const SGResidualOp&);

  protected:

    //! Number of blocks
    unsigned int nblock;

    //! Map for base discretization
    Teuchos::RCP<const Epetra_Map> map;

    //! Stochastic Galerking basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > sg_basis;

    //! SG time derivative vector (may be null)
    Teuchos::RCP<const EpetraExt::BlockVector> sg_xdot;

    //! SG solution vector
    Teuchos::RCP<const EpetraExt::BlockVector> sg_x;

    //! SG residual vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_f;

    //! Time derivative vector (may be null)
    std::vector< Teuchos::RCP<Epetra_Vector> > xdot;

    //! Solution vector
    std::vector< Teuchos::RCP<Epetra_Vector> > x;

  };

  //! Fill operator for Jacobian
  class SGJacobianOp : 
    public FEApp::AbstractInitPostOp<FEApp::SGJacobianType> {
  public:

    typedef SGType::expansion_type::tp_type tp_type;

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    SGJacobianOp(
         const Teuchos::RCP<const Epetra_Map>& base_map,
         const Teuchos::RCP<const Epetra_CrsGraph>& base_graph,
	       const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis,
         const Teuchos::RCP<const tp_type >& Cijk,
         double alpha, double beta,
         const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
         const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
         const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f,
         const Teuchos::RCP<EpetraExt::BlockCrsMatrix>& sg_overlapped_jac);

    //! Destructor
    virtual ~SGJacobianOp();

    //! Reset operator for new fill
    virtual void reset(
          double alpha, double beta,
	        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
          const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x);

    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< SGFadType >* elem_xdot,
                             std::vector< SGFadType >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< SGFadType >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< SGFadType >* node_xdot,
                          std::vector< SGFadType >& node_x);
    
    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< SGFadType >& node_f);

    //! Finalize fill
    virtual void finalizeFill();

  private:
    
    //! Private to prohibit copying
    SGJacobianOp(const SGJacobianOp&);
    
    //! Private to prohibit copying
    SGJacobianOp& operator=(const SGJacobianOp&);

  protected:

    //! Number of blocks
    unsigned int nblock;

    //! Map for base discretization
    Teuchos::RCP<const Epetra_Map> map;

    //! Graph for base discretization
    Teuchos::RCP<const Epetra_CrsGraph> graph;

    //! Stochastic Galerkin basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > sg_basis;

    //! Stochastic Galerkin triple product
    Teuchos::RCP<const tp_type> Cijk;

    //! Coefficient of mass matrix
    double m_coeff;

    //! Coefficient of Jacobian matrix
    double j_coeff;

    //! SG time derivative vector (may be null)
    Teuchos::RCP<const EpetraExt::BlockVector> sg_xdot;

    //! SG solution vector
    Teuchos::RCP<const EpetraExt::BlockVector> sg_x;

    //! SG residual vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_f;

    //! SG Jacobian matrix
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> sg_jac;

    //! Time derivative vector (may be null)
    std::vector< Teuchos::RCP<Epetra_Vector> > xdot;

    //! Solution vector
    std::vector< Teuchos::RCP<Epetra_Vector> > x;

  };

    //! Fill operator for Jacobian
  class SGMatrixFreeJacobianOp : 
    public FEApp::AbstractInitPostOp<FEApp::SGJacobianType> {
  public:

    typedef SGType::expansion_type::tp_type tp_type;

    //! Constructor
    /*!
     * Set xdot to Teuchos::null for steady-state problems
     */
    SGMatrixFreeJacobianOp(
         const Teuchos::RCP<const Epetra_Map>& base_map,
         const Teuchos::RCP<const Epetra_CrsGraph>& base_graph,
	       const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis,
         double alpha, double beta,
         const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
         const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
         const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f);

    //! Destructor
    virtual ~SGMatrixFreeJacobianOp();

    //! Reset operator for new fill
    virtual void reset(
	        double alpha, double beta,
	        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
          const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x);

    //! Get Jacobian blocks
    virtual std::vector< Teuchos::RCP<Epetra_CrsMatrix> >&
    getJacobianBlocks();

    //! Evaulate element init operator
    virtual void elementInit(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< SGFadType >* elem_xdot,
                             std::vector< SGFadType >& elem_x);

    //! Evaluate element post operator
    virtual void elementPost(const FEApp::AbstractElement& e,
                             unsigned int neqn,
                             std::vector< SGFadType >& elem_f);

    //! Evaulate node init operator
    virtual void nodeInit(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< SGFadType >* node_xdot,
                          std::vector< SGFadType >& node_x);

    //! Evaluate node post operator
    virtual void nodePost(const FEApp::NodeBC& bc,
                          unsigned int neqn,
                          std::vector< SGFadType >& node_f);

    //! Finalize fill
    virtual void finalizeFill();

  private:
    
    //! Private to prohibit copying
    SGMatrixFreeJacobianOp(const SGMatrixFreeJacobianOp&);
    
    //! Private to prohibit copying
    SGMatrixFreeJacobianOp& operator=(const SGMatrixFreeJacobianOp&);

  protected:

    //! Number of blocks
    unsigned int nblock;

    //! Map for base discretization
    Teuchos::RCP<const Epetra_Map> map;

    //! Graph for base discretization
    Teuchos::RCP<const Epetra_CrsGraph> graph;

    //! Stochastic Galerkin basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> > sg_basis;

    //! Coefficient of mass matrix
    double m_coeff;

    //! Coefficient of Jacobian matrix
    double j_coeff;

    //! SG time derivative vector (may be null)
    Teuchos::RCP<const EpetraExt::BlockVector> sg_xdot;

    //! SG solution vector
    Teuchos::RCP<const EpetraExt::BlockVector> sg_x;

    //! SG residual vector
    Teuchos::RCP<EpetraExt::BlockVector> sg_f;

    //! SG Jacobian matrix
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> sg_jac;

    //! Time derivative vector (may be null)
    std::vector< Teuchos::RCP<Epetra_Vector> > xdot;

    //! Solution vector
    std::vector< Teuchos::RCP<Epetra_Vector> > x;

    //! Jacobian matrix
    std::vector< Teuchos::RCP<Epetra_CrsMatrix> > jac;

  };

#endif // SG_ACTIVE

}

#endif // FEAPP_INITPOSTOPS_HPP
