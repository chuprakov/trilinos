/*
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
*/

#ifndef EPETRA2TSFUTILS_H
#define EPETRA2TSFUTILS_H

#include "kay.h"

#define EPETRA_INVERSE 1
#define EPETRA_MATRIX  2

using namespace TSF;

  TSFVectorSpace EpetraCRS2TSFVspace(Epetra_CrsMatrix *Mx_crs);
  int* Vbr2TSF(int Nvel, int Npress);
  TSFLinearOperator EpetraCRS2TSF(TSFVectorSpace vSpace,TSFVectorSpace pSpace,Epetra_CrsMatrix *Mx_crs);

#endif


