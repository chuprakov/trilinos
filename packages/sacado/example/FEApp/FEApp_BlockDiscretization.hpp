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

#ifndef FEAPP_BLOCKDISCRETIZATION_HPP
#define FEAPP_BLOCKDISCRETIZATION_HPP

#include <vector>
#include "FEApp_TemplateTypes.hpp"

#if SG_ACTIVE

#include "Epetra_Comm.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "EpetraExt_BlockUtility.h"

#include "FEApp_AbstractDiscretization.hpp"

#include "Stokhos_OrthogPolyBasis.hpp"

namespace FEApp {

  class BlockDiscretization : public FEApp::AbstractDiscretization {
  public:

    //! Constructor
    BlockDiscretization(
     const Teuchos::RCP<const Epetra_Comm>& comm,
     const Teuchos::RCP<const FEApp::AbstractDiscretization>& underlyingDisc_,
     const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
     bool makeJacobian);

    //! Destructor
    virtual ~BlockDiscretization();

    //! Create element mesh
    virtual void createMesh();

    //! Create DOF maps
    virtual void createMaps();

    //! Create Jacobian graph
    virtual void createJacobianGraphs();

    //! Get element mesh
    virtual Teuchos::RCP<const FEApp::Mesh> 
    getMesh() const; 

    //! Get DOF map
    virtual Teuchos::RCP<const Epetra_Map> 
    getMap() const;

    //! Get overlapped DOF map
    virtual Teuchos::RCP<const Epetra_Map> 
    getOverlapMap() const;

    //! Get Jacobian graph
    virtual Teuchos::RCP<const Epetra_CrsGraph> 
    getJacobianGraph() const;

    //! Get overlap Jacobian graph
    virtual Teuchos::RCP<const Epetra_CrsGraph> 
    getOverlapJacobianGraph() const;

    //! Get number of nodes per element
    virtual int getNumNodesPerElement() const;

    //! Get Jacobian
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> 
    getJacobian();

    //! Get overlap Jacobian
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> 
    getOverlapJacobian();


  private:

    //! Private to prohibit copying
    BlockDiscretization(const BlockDiscretization&);

    //! Private to prohibit copying
    BlockDiscretization& operator=(const BlockDiscretization&);

  protected:
    
    //! Underlying discretization object
    Teuchos::RCP<const FEApp::AbstractDiscretization> underlyingDisc;

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> globalComm;

    //! Stochastic Galerkin basis
    Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > sg_basis;

    //! Unknown Map
    Teuchos::RCP<const Epetra_Map> map;

    //! Overlapped unknown map
    Teuchos::RCP<const Epetra_Map> overlap_map;

    //! Jacobian matrix graph
    Teuchos::RCP<const Epetra_CrsGraph> graph;

    //! Overlapped Jacobian matrix graph
    Teuchos::RCP<const Epetra_CrsGraph> overlap_graph;

    //! Jacobian matrix
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> jac;

    //! Overlapped Jacobian matrix
    Teuchos::RCP<EpetraExt::BlockCrsMatrix> overlap_jac;

  };

}

#endif // SG_ACTIVE

#endif // FEAPP_CZERODISCRETIZATION_HPP
