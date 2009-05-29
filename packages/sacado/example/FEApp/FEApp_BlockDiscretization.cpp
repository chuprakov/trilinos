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
#include "Epetra_Export.h"

#include "FEApp_BlockDiscretization.hpp"

#ifdef HAVE_MPI
#include "EpetraExt_MultiMpiComm.h"
#else
#include "EpetraExt_MultiSerialComm.h"
#endif

#if SG_ACTIVE

FEApp::BlockDiscretization::BlockDiscretization(
    const Teuchos::RCP<const Epetra_Comm>& comm,
    const Teuchos::RCP<const FEApp::AbstractDiscretization>& underlyingDisc_,
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis_,
     bool makeJacobian) :
  underlyingDisc(underlyingDisc_),
  sg_basis(sg_basis_)
{

  unsigned int num_sg_blocks = sg_basis->size();
#ifdef HAVE_MPI
  // No parallelism over blocks, so spatial partition is unchanged 
  // as comm->NumProc()
  Teuchos::RCP<EpetraExt::MultiComm> multiComm =
    Teuchos::rcp(new EpetraExt::MultiMpiComm(MPI_COMM_WORLD, 
					     comm->NumProc(), 
					     num_sg_blocks));
#else
  Teuchos::RCP<EpetraExt::MultiComm> multiComm =
    Teuchos::rcp(new EpetraExt::MultiSerialComm(num_sg_blocks));
#endif

  // Create block matrix and graph from underlyingDisc and
  // the block information stored in globalComm
  int numBlockRows =  multiComm->NumTimeSteps();
  int myBlockRows  =  multiComm->NumTimeStepsOnDomain();
  int myFirstBlockRow = multiComm->FirstTimeStepOnDomain();
  globalComm = multiComm;

  // DENSE STENCIL for Stochastic Galerkin
  // For 3 blocks on 2 procs, this should be:
  // Proc  nBR  mBR  mFBR     Stencil      Index
  //  0     3    2    0       0  1  2        0
  //                         -1  0  1        1
  //  1     3    1    2      -2 -1  0        2
  //
   std::vector< std::vector<int> > rowStencil(myBlockRows);
   std::vector<int> rowIndex(myBlockRows);
   for (int i=0; i < myBlockRows; i++) {
     for (int j=0; j < numBlockRows; j++) 
       rowStencil[i].push_back(-myFirstBlockRow - i + j);
     rowIndex[i] = (i + myFirstBlockRow);
   }

   if (makeJacobian) {
     //Construct BlockGraphs from underlying graphs and stencils
     jac = Teuchos::rcp(new EpetraExt::BlockCrsMatrix(
				  *(underlyingDisc->getJacobianGraph()),
				  rowStencil, rowIndex, *globalComm));
     overlap_jac = Teuchos::rcp(new EpetraExt::BlockCrsMatrix(
				  *(underlyingDisc->getOverlapJacobianGraph()),
				  rowStencil, rowIndex, *globalComm));
     graph = Teuchos::rcp(&(jac->Graph()),false);
     overlap_graph = Teuchos::rcp(&(overlap_jac->Graph()),false);
     map = Teuchos::rcp(&(jac->RowMatrixRowMap()),false);
     overlap_map = Teuchos::rcp(&(overlap_jac->RowMatrixRowMap()),false);
   }
   else {
     map = Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
						  *(underlyingDisc->getMap()),
						  rowIndex,
						  *globalComm));
     overlap_map = Teuchos::rcp(EpetraExt::BlockUtility::GenerateBlockMap(
					    *(underlyingDisc->getOverlapMap()),
					    rowIndex,
					    *globalComm));
   }
}

FEApp::BlockDiscretization::~BlockDiscretization()
{
}

void
FEApp::BlockDiscretization::createMesh()
{
  // Done in underlying Discretization
}

void
FEApp::BlockDiscretization::createMaps()
{
  // Done in constructor
}

void
FEApp::BlockDiscretization::createJacobianGraphs()
{
  // Done in constructor
}
	    
Teuchos::RCP<const FEApp::Mesh>
FEApp::BlockDiscretization::getMesh() const
{
  return underlyingDisc->getMesh();
}

Teuchos::RCP<const Epetra_Map>
FEApp::BlockDiscretization::getMap() const
{
  return map;
}

Teuchos::RCP<const Epetra_Map>
FEApp::BlockDiscretization::getOverlapMap() const
{
  return overlap_map;
}

Teuchos::RCP<const Epetra_CrsGraph>
FEApp::BlockDiscretization::getJacobianGraph() const
{
  return graph;
}

Teuchos::RCP<const Epetra_CrsGraph>
FEApp::BlockDiscretization::getOverlapJacobianGraph() const
{
  return overlap_graph;
}

int
FEApp::BlockDiscretization::getNumNodesPerElement() const
{
  return underlyingDisc->getNumNodesPerElement();
}

Teuchos::RCP<EpetraExt::BlockCrsMatrix> 
FEApp::BlockDiscretization::getJacobian()
{
  return jac;
}

Teuchos::RCP<EpetraExt::BlockCrsMatrix> 
FEApp::BlockDiscretization::getOverlapJacobian()
{
  return overlap_jac;
}

#endif // SG_ACTIVE
