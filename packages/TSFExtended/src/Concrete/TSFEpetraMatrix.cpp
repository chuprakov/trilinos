/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
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
// **********************************************************************/
/* @HEADER@ */

#include "TSFEpetraMatrix.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_MPIComm.hpp"

using namespace TSFExtended;
using namespace Teuchos;

EpetraMatrix::EpetraMatrix(const RefCountPtr<const EpetraVectorSpace>& domain,
                           const RefCountPtr<const EpetraVectorSpace>& range)
  : TSFCore::EpetraLinearOp()
{
  /* initializing ncols to zero allows later fill */
  RefCountPtr<Epetra_CrsMatrix> A 
    = rcp(new Epetra_CrsMatrix(Copy, *(range->epetra_map()), 0));

  initialize(A, TSFCore::NOTRANS);
}

void EpetraMatrix::setGraph(int nLocalRows,
                            const int* globalRowIndex,
                            const int* numNonzeros,
                            const int** columnIndices)
{
  Array<double> zeros;
  int z = zeros.size();

  for (int i=0; i<nLocalRows; i++)
    {
      int g = globalRowIndex[i];
      int nnz = numNonzeros[i];
      const int* col = columnIndices[i];
      if (nnz > z)
        {
          zeros.resize(nnz);
          for (int j=z; j<nnz; j++) zeros[i] = 0.0;
        }
      int ierr = crsMatrix()->InsertGlobalValues(g, nnz, 
                                                 &(zeros[0]), 
                                                 (int*) col);
      TEST_FOR_EXCEPTION(ierr < 0, runtime_error,
                         "failed to configure row " << g 
                         << " in EpetraMatrix::setGraph() with nnz="
                         << nnz << ". Error code was " << ierr);
    }
}

void EpetraMatrix::freezeValues()
{
  int ierr = crsMatrix()->FillComplete();

  TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                     "EpetraMatrix::freezeValues() failed during call "
                     "to FillComplete(). Error code was " << ierr);
}

void EpetraMatrix::addToRow(int globalRowIndex,
                            int nElemsToInsert,
                            const int* globalColumnIndices,
                            const double* elementValues)
{
  int ierr = crsMatrix()->SumIntoGlobalValues(globalRowIndex,
                                              nElemsToInsert,
                                              (double*) elementValues,
                                              (int*) globalColumnIndices);

  TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                     "failed to add to row " << globalRowIndex
                     << " in EpetraMatrix::addToRow() with nnz="
                     << nElemsToInsert 
                     << ". Error code was " << ierr);
}

void EpetraMatrix::setRowValues(int globalRowIndex,
                                int nElemsToInsert,
                                const int* globalColumnIndices,
                                const double* elementValues)
{
  int ierr = crsMatrix()->InsertGlobalValues(globalRowIndex,
                                             nElemsToInsert,
                                             (double*) elementValues,
                                             (int*) globalColumnIndices);

  TEST_FOR_EXCEPTION(ierr < 0, runtime_error, 
                     "failed to add to row " << globalRowIndex
                     << " in EpetraMatrix::setRowValues() with nnz="
                     << nElemsToInsert 
                     << ". Error code was " << ierr);
}

void EpetraMatrix::zero()
{
  crsMatrix()->PutScalar(0.0);
}

void EpetraMatrix::print(ostream& os) const 
{
  int nProc = MPISession::getNProc();
  int rank = MPISession::getRank();
  for (int i=0; i<nProc; i++)
    {
      MPIComm::world().synchronize();
      if (i==rank) crsMatrix()->Print(os);
    }
}

string EpetraMatrix::describe() const 
{
  string rtn = "EpetraMatrix[nRow=" 
    + Teuchos::toString(crsMatrix()->NumGlobalRows())
    + ", nCol=" + Teuchos::toString(crsMatrix()->NumGlobalCols())
    + "]";
  return rtn;
}

Epetra_CrsMatrix* EpetraMatrix::crsMatrix()
{
  Epetra_CrsMatrix* crs 
    = dynamic_cast<Epetra_CrsMatrix*>(epetra_op().get());

  TEST_FOR_EXCEPTION(crs==0, runtime_error,
                     "cast failed in EpetraMatrix::crsMatrix()");

  return crs;
}

const Epetra_CrsMatrix* EpetraMatrix::crsMatrix() const 
{
  const Epetra_CrsMatrix* crs 
    = dynamic_cast<const Epetra_CrsMatrix*>(epetra_op().get());

  TEST_FOR_EXCEPTION(crs==0, runtime_error,
                     "cast failed in EpetraMatrix::crsMatrix()");

  return crs;
}

RefCountPtr<const TSFCore::EpetraVectorSpace> 
EpetraMatrix::allocateDomain(const RefCountPtr<Epetra_Operator>  &op 
                               ,TSFCore::ETransp  op_trans 
                               )  const
{
  return rcp( new EpetraVectorSpace( rcp(&op->OperatorDomainMap(),false) ) );
}

RefCountPtr<const TSFCore::EpetraVectorSpace> 
EpetraMatrix::allocateRange(const RefCountPtr<Epetra_Operator>  &op 
                              ,TSFCore::ETransp  op_trans 
                               )  const
{
  return rcp( new EpetraVectorSpace( rcp(&op->OperatorRangeMap(),false) ) );
}

Epetra_CrsMatrix& EpetraMatrix::getConcrete(const LinearOperator<double>& A)
{
  EpetraMatrix* ep 
    = dynamic_cast<EpetraMatrix*>(A.ptr().get());
  TEST_FOR_EXCEPTION(ep==0, std::runtime_error,
                     "EpetraMatrix::getConcrete called on a matrix that "
                     "could not be cast to an EpetraMatrix");
  return *(ep->crsMatrix());
}
