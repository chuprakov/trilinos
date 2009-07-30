//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER

#include "euclid_Helpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include <iostream>
#include <fstream>
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"
#include "Ifpack_Hypre.h"
#include "Epetra_MpiComm.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <sstream>
#include <map>

using Teuchos::RCP;
using Teuchos::rcp;

Epetra_CrsMatrix::Epetra_CrsMatrix* newCrsMatrix(int N){
  int ierr = 0;
  Epetra_MpiComm Comm(MPI_COMM_WORLD);

  RCP<Epetra_Map>            Map;
  // pointer to the matrix to be created
  Epetra_CrsMatrix*      Matrix;
  // container for parameters
  /*Teuchos::ParameterList GaleriList;
  int nx = N * Comm.NumProc();
  int ny = N * Comm.NumProc();
  GaleriList.set("nx", nx);
  GaleriList.set("ny", ny);

  Map = rcp(Galeri::CreateMap("Cartesian2D", Comm, GaleriList));
  Matrix   = Galeri::CreateCrsMatrix("Laplace2D", Map.get(), GaleriList);*/
  
  Map = rcp(new Epetra_Map(N*Comm.NumProc(), N, 0, Comm));
  Matrix = new Epetra_CrsMatrix(Copy, *Map, 3);

  int ilower = N*Comm.MyPID(); int iupper = N*(Comm.MyPID()+1)-1;
  for(int i = ilower; i <= iupper; i++){
    if(i != 0 && i != N*Comm.NumProc()-1){
      Teuchos::Array<int> indices; 
      Teuchos::Array<double> values;
      values.resize(3); indices.resize(3);
      indices[0] = i-1;
      indices[1] = i;
      indices[2] = i+1;
      values[0] = -1.0;
      values[1] = 4.0;
      values[2] = -1.0;
      Matrix->InsertGlobalValues(i, 3, &values[0], &indices[0]);
    } else if(i != N*Comm.NumProc()-1){
      Teuchos::Array<int> indices;
      Teuchos::Array<double> values;
      values.resize(2); indices.resize(2);
      indices[0] = 0;
      indices[1] = 1;
      values[0] = 4.0;
      values[1] = -1.0;
      Matrix->InsertGlobalValues(0, 2, &values[0], &indices[0]);
    } else{
      Teuchos::Array<int> indices;
      Teuchos::Array<double> values;
      values.resize(2); indices.resize(2);
      indices[0] = i-1;
      indices[1] = i;
      values[0] = -1.0;
      values[1] = 4.0;
      Matrix->InsertGlobalValues(i, 2, &values[0], &indices[0]);
    }
  }

  ierr += Matrix->FillComplete();
  if(ierr != 0){
    printf("\n\nMatrix was not created properly!\n\n");
  }
  return Matrix;
}


bool EquivalentVectors(Epetra_MultiVector &Y1, Epetra_MultiVector &Y2, const double tol){
  
  bool retVal = true;
  
  int num_vectors = Y1.NumVectors();
  if(Y2.NumVectors() != num_vectors){
    printf("Multivectors do not have same number of vectors.\n");
    return false;
  }
  
  for(int j = 0; j < num_vectors; j++){
    if(Y1.MyLength() != Y2.MyLength()){
      printf("Vectors are not same size on local processor.\n");
      return false;
    }
    for(int i = 0; i < Y1.GlobalLength(); i++){
      int Y1_local = Y1.Map().LID(i);
      int Y2_local = Y2.Map().LID(i);
      if(Y1_local < 0 || Y2_local < 0){
        continue;
      }
      if(fabs((*Y1(j))[Y1_local] - (*Y2(j))[Y2_local]) > tol){
        printf("Vector number[%d] ", j);
        printf("Val1[%d] = %f != Val2[%d] = %f\n", i, (*Y1(j))[Y1_local], i, (*Y2(j))[Y2_local]);
        retVal = false;
      }
    }
  }
  Teuchos::Array<int> vals; vals.resize(Y1.Comm().NumProc());
  int my_vals[1]; my_vals[0] = (int)retVal;
  Y1.Comm().GatherAll(my_vals, &vals[0], 1);
  for(int i = 0; i < Y1.Comm().NumProc(); i++){
    if(vals[i] == false){
      retVal = false;
    }
  }
  if(retVal == false){
    printf("[%d]Failed vector equivalency test.\n", Y1.Comm().MyPID());
  }
  return retVal;
}

