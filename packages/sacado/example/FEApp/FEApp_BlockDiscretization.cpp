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

#if SG_ACTIVE

FEApp::BlockDiscretization::BlockDiscretization(
    const Teuchos::RCP<const EpetraExt::MultiMpiComm>& globalComm_,
    const Teuchos::RCP<const FEApp::AbstractDiscretization>& underlyingDisc_,
    const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis_ ) :
  underlyingDisc(underlyingDisc_),
  globalComm(globalComm_),
  sg_basis(sg_basis_)
{

  // Create block matrix and graph from underlyingDisc and
  // the block information stored in globalComm

  int numBlockRows =  globalComm->NumTimeSteps();
  int myBlockRows  =  globalComm->NumTimeStepsOnDomain();
  int myFirstBlockRow = globalComm->FirstTimeStepOnDomain();

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

//    graph = Teuchos::RCP<Epetra_CrsGraph>(
//        EpetraExt::BlockUtility::GenerateBlockGraph(
//                   *(underlyingDisc->getJacobianGraph()),
//                   rowStencil, rowIndex, *globalComm));

//    overlap_graph = Teuchos::RCP<Epetra_CrsGraph>(
//        EpetraExt::BlockUtility::GenerateBlockGraph(
//                   *(underlyingDisc->getOverlapJacobianGraph()),
//                   rowStencil, rowIndex, *globalComm));

//    //Temporarily constructing CrsMatrix from graph, so the Epetra_Map
//    // can be pulled out. I only see how to get Epetra_BlockMap from the graph

//    Epetra_CrsMatrix* mat;

//    mat = Teuchos::rcp(new EpetraExt::BlockCrsMatrix(View, *graph);
//    map =  Teuchos::RCP<Epetra_Map>(new Epetra_Map(mat->RowMatrixRowMap()));
//    delete mat;

//    mat = new Epetra_CrsMatrix(View, *overlap_graph);
//    overlap_map =  Teuchos::RCP<Epetra_Map>(new Epetra_Map(mat->RowMatrixRowMap()));
//    delete mat;

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
