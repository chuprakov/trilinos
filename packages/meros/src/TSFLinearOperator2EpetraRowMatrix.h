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

#ifndef TSFLINEAROPERATOR2EPETRAROWMATRIX_H
#define TSFLINEAROPERATOR2EPETRAROWMATRIX_H

#include "Epetra_SerialComm.h"
#include "Epetra_Comm.h"
#include "Epetra_Vector.h"
#include "TSFMPI.h"
#include "TSFOut.h"
#include "TSFVectorType.h"
#include "TSFProductSpace.h"
#include "TSFLinearOperator.h"
#include "TSFBlockLinearOperator.h"
#include "PetraVectorType.h"
#include "PetraVectorSpace.h"
#include "PetraVector.h"
#include "PetraMatrix.h"
#include "GMRESSolver.h"
#include "AZTECSolver.h"
#include "TSFPreconditionerFactory.h"
#include "GenericRightPreconditioner.h"
#include "TSFPreconditioner.h"
#include "TSFMatrixOperator.h"

#define EPETRA_INVERSE 1
#define EPETRA_MATRIX  2

// Forward declaration
class Epetra_Vector;
using namespace TSF;

class TSFLinearOperator2EpetraRowMatrix : public virtual Epetra_RowMatrix  {

 public:

  TSFLinearOperator2EpetraRowMatrix(TSFLinearOperator A_tsf, const Epetra_Comm *Comm_pet,
                                    Epetra_Map *Map_pet, int *map, int matrix_type);

  TSFLinearOperator2EpetraRowMatrix(const Epetra_Comm *Comm_pet, int matrix_type);
  ~TSFLinearOperator2EpetraRowMatrix();

  int NumMyRows() const;
  int NumMyCols() const;
  double NormInf() const;
  bool HasNormInf() const;
  const Epetra_Map & OperatorDomainMap() const;
  const Epetra_Map & OperatorRangeMap() const;
  const Epetra_Comm & Comm() const;
  int SetUseTranspose(bool UseTranspose);
  char  *Label() const;

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //***************************************************************************
  // NOT CALLED BY US!!!!
  //***************************************************************************
  int MaxNumEntries () const;

  const Epetra_BlockMap &Map () const;

  int NumMyRowEntries(int MyRow, int & NumEntries) const;

  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

  int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const;

  int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X,
            Epetra_MultiVector& Y) const;

  int InvRowSums(Epetra_Vector& x) const;

  int LeftScale(const Epetra_Vector& x);

  int InvColSums(Epetra_Vector& x) const;

  bool    Filled() const;
  double NormOne() const;
  int NumGlobalNonzeros() const;
  int NumGlobalRows() const;
  int NumGlobalCols() const;
  int NumGlobalDiagonals() const;
  int NumMyNonzeros() const;
  int NumMyDiagonals() const;
  bool LowerTriangular() const;
  bool UpperTriangular() const;
  bool UseTranspose() const;

  int RightScale(const Epetra_Vector& x);

  const Epetra_Map & RowMatrixRowMap() const;

  const Epetra_Map & RowMatrixColMap() const;

  const Epetra_Import * RowMatrixImporter() const;


  // The current applyoperator() is very specific to the needs of MPSalsa. In particular,
  // it is assumed that the TSF matrix that is being wrapped is a block matrix that
  // operates on a series of petra type vector spaces. It is also assumed that the
  // data needs to be reshuffled in the vectors. For example, the 'X' and 'Y' vectors
  // correspond to something like
  //                ( u0, v0, p0, u1, v1, p1, u2, v2, p2, .... )
  // and we want to reorder this like
  //                ( u0, v0, u1,  v1, u2, v2, .... )  and (p0, p1, p2, ...)


  int applyoperator(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


  TSFLinearOperator getTSF( );
  int *getBlockAssignments();


  int initialize(TSFLinearOperator A_tsf,  Epetra_Map *Map_pet, int *map);

 protected:


  TSFVectorSpace rangeBlockSpace_, domainBlockSpace_;
  int numRangeBlocks_, numDomainBlocks_, rangeLength_, domainLength_;
  int matrix_type_;

  int *map_;
  bool OwnMap_;
  TSFLinearOperator Matrix_;
  const Epetra_Map *OperatorDomainMap_;
  const Epetra_Comm *Comm_;
  Epetra_Import *dummy;
  char  *myname;
  char  *nop;
};
#endif
