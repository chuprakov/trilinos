//@HEADER
// ************************************************************************
// 
//               Tpetra: Templated Linear Algebra Services Package 
//                 Copyright (2001) Sandia Corporation
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
// ************************************************************************
//@HEADER

#include "Tpetra_ConfigDefs.hpp"
#ifdef HAVE_MPI
#include "Tpetra_MpiPlatform.hpp"
#include "Tpetra_MpiComm.hpp"
#else
#include "Tpetra_SerialPlatform.hpp"
#include "Tpetra_SerialComm.hpp"
#endif
#include "Tpetra_ElementSpace.hpp"
#include "Tpetra_VectorSpace.hpp"
#include "Tpetra_CisMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Import.hpp"
#include "Teuchos_ScalarTraits.hpp"

typedef int OrdinalType;
typedef double ScalarType;

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Tpetra::MpiComm<OrdinalType, ScalarType> Comm(MPI_COMM_WORLD);
#else
  Tpetra::SerialComm<OrdinalType, ScalarType> Comm;
#endif

  if (Comm.getNumImages() != 2)
  {
    cout << "This example must be ran with two processors" << endl;
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
#endif
  }

  // Get zero and one for the OrdinalType
  
  OrdinalType const OrdinalZero = Teuchos::ScalarTraits<OrdinalType>::zero();
  OrdinalType const OrdinalOne  = Teuchos::ScalarTraits<OrdinalType>::one();

  // Get zero and one for the ScalarType
  
  ScalarType const ScalarZero = Teuchos::ScalarTraits<ScalarType>::zero();
  ScalarType const ScalarOne  = Teuchos::ScalarTraits<ScalarType>::one();

  // Creates a vector of size `length', then set the elements values.
  
  OrdinalType length    = OrdinalOne * 10;
  OrdinalType indexBase = OrdinalZero;

  // Creation of a platform
  
#ifdef HAVE_MPI
  const Tpetra::MpiPlatform <OrdinalType, OrdinalType> platformE(MPI_COMM_WORLD);
  const Tpetra::MpiPlatform <OrdinalType, ScalarType> platformV(MPI_COMM_WORLD);
#else
  const Tpetra::SerialPlatform <OrdinalType, OrdinalType> platformE;
  const Tpetra::SerialPlatform <OrdinalType, ScalarType> platformV;
#endif

  // These are the `elements' assigned to this image
  
  OrdinalType NumMyElements = OrdinalOne * 3;
  vector<OrdinalType> MyGlobalElements(NumMyElements);

  if (Comm.getMyImageID() == 0)
  {
    for (OrdinalType i = OrdinalZero ; i < NumMyElements ; ++i)
      MyGlobalElements[OrdinalZero + i] = OrdinalZero + i;
  }
  else
  {
    for (OrdinalType i = OrdinalZero ; i < NumMyElements ; ++i)
      MyGlobalElements[OrdinalZero + i] = OrdinalZero + i + 3;
  }

  Tpetra::ElementSpace<OrdinalType> 
    elementSpace(-1, NumMyElements, MyGlobalElements, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> vectorSpace(elementSpace, platformV);

  cout << elementSpace;

  // These are the `elements' + ghost nodes

  OrdinalType NumMyPaddedElements = OrdinalOne * 4;
  vector<OrdinalType> MyGlobalPaddedElements(NumMyPaddedElements);

  if (Comm.getMyImageID() == 0)
  {
    MyGlobalPaddedElements[0] = 0;
    MyGlobalPaddedElements[1] = 1;
    MyGlobalPaddedElements[2] = 2;
    MyGlobalPaddedElements[3] = 3;
  }
  else
  {
    MyGlobalPaddedElements[0] = 3;
    MyGlobalPaddedElements[1] = 4;
    MyGlobalPaddedElements[2] = 5;
    MyGlobalPaddedElements[3] = 2;
  }

  Tpetra::ElementSpace<OrdinalType> 
    elementPaddedSpace(-1, NumMyPaddedElements, MyGlobalPaddedElements, indexBase, platformE);
  Tpetra::VectorSpace<OrdinalType, ScalarType> vectorPaddedSpace(elementPaddedSpace, platformV);

  cout << elementPaddedSpace;

  // This is the connectivity, in local numbering

  OrdinalType NumMyFiniteElements = OrdinalOne * 3;
  OrdinalType Connectivity[NumMyFiniteElements][2];
  ScalarType  Coord[NumMyPaddedElements][2];

  if (Comm.getMyImageID() == 0)
  {
    //     e0     e1     e2
    // (0) -- (1) -- (2) -- [3] -- (X) -- (X)
    Connectivity[0][0] = 0; Connectivity[0][1] = 1;
    Connectivity[1][0] = 1; Connectivity[1][1] = 2;
    Connectivity[2][0] = 2; Connectivity[2][1] = 3;
  }
  else
  {
    //                   e2     e0     e1
    // (X) -- (X) -- [3] -- (0) -- (1) -- (2)
    Connectivity[0][0] = 0; Connectivity[0][1] = 1;
    Connectivity[1][0] = 1; Connectivity[1][1] = 2;
    Connectivity[2][0] = 3; Connectivity[2][1] = 0;
  }

  // Setup the matrix
  
  Tpetra::CisMatrix<OrdinalType,ScalarType> matrix(vectorSpace);

  for (OrdinalType FEID = OrdinalZero ; FEID < NumMyFiniteElements ; ++FEID)
  {
    cout << "element = " << FEID << endl;
    // build the local matrix, in this case A_loc
    ScalarType A_loc[2][2] = {ScalarOne, -ScalarOne, -ScalarOne, ScalarOne};

    vector<OrdinalType> GIDs(2);

    // submit entries of A_loc into matrix
    for (OrdinalType LRID = OrdinalZero ; LRID < 2 ; ++LRID)
    {
      cout << "-(" << Connectivity[FEID][LRID] << ")-";
      GIDs[LRID] = elementPaddedSpace.getGID(Connectivity[FEID][LRID]);
      cout << GIDs[LRID] << endl;
    }

    for (OrdinalType LRID = OrdinalZero ; LRID < 2 ; ++LRID)
    {
      if (Connectivity[FEID][LRID] < NumMyElements)
      {
        for (OrdinalType LCID = OrdinalZero ; LCID < 2 ; ++LCID)
        {
          if (GIDs[LRID] != -OrdinalOne)
            matrix.submitEntry(Tpetra::Add, GIDs[LRID], A_loc[LRID][LCID],
                               GIDs[LCID]);
        }
      }
    }
  }

  matrix.fillComplete();

  cout << matrix;
#ifdef HAVE_MPI
  MPI_Finalize() ;
#endif

  return(EXIT_SUCCESS);
}
