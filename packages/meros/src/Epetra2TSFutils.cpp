// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Epetra2TSFutils.h"

// Forward declaration
using namespace TSF;

// Build a map. This map is intended to go between vbr-style vectors and TSF-style vectors. 
 // Since this example does not contain VBR data, the map is pointless but needed
 // for the proper creation of an epetra wrapper around a TSF matrix.
 // Specifically, map[i] != 1 indicates that the ith variable is in the
 //                           first block.
 //               map[i] == 1 indicates that the ith variable is in the
 //                           second block.
int* Vbr2TSF(int Nvel, int Npress)
{
 int *tmp = (int *) malloc(sizeof(int)*(Nvel + Npress));
 for (int i = 0; i < Nvel; i++) tmp[i] = 0;
 for (int i = Nvel; i < Nvel + Npress; i++)  tmp[i] = 1;
 return tmp;
}

// Convert block matrices to TSF matrices.
TSFVectorSpace EpetraCRS2TSFVspace(Epetra_CrsMatrix *Mx_crs)
{
  Epetra_Map *p1 = (Epetra_Map *) &(Mx_crs->OperatorDomainMap());
  Epetra_Map *p2 = (Epetra_Map *) &(Mx_crs->RowMatrixColMap());
  Epetra_Import *p3 = (Epetra_Import *) Mx_crs->RowMatrixImporter();
  TSFSmartPtr<Epetra_Map> tp1 = TSFSmartPtr<Epetra_Map>(p1,false);
  TSFSmartPtr<Epetra_Map> tp2 = TSFSmartPtr<Epetra_Map>(p2,false);
  TSFSmartPtr<Epetra_Import> tp3 = TSFSmartPtr<Epetra_Import>(p3,false);
  TSFVectorSpace tmpSpace = new PetraVectorSpace( tp1, tp2, tp3);
  return tmpSpace;
}

// Build TSF vector spaces
// Be sure to call with TSF style space ordering (i.e., domain, range)
TSFLinearOperator EpetraCRS2TSF(TSFVectorSpace vSpace,TSFVectorSpace pSpace,Epetra_CrsMatrix *Mx_crs)
{
 PetraMatrix* Mx_petra = new PetraMatrix(vSpace, pSpace);
 Mx_petra->setPetraMatrix(Mx_crs,true); // insert epetra matrix into TSF matrix.
 TSFLinearOperator Mx_tsf = Mx_petra;   // 'true' indicates that when this TSF matrix
                                         // is deleted (automatically via smart
                                         // pointers),  the underlying epetra matrix is 
                                         // also deleted.
 return Mx_tsf;
}

 
