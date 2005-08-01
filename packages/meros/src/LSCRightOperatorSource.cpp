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

#include "LSCRightOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFBlockLinearOperator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "PetraMatrix.h"

using namespace TSF;
using namespace Meros;
using std::string;


LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S)
  : S_(S), Ap_(), hasAp_(false), Dinv_(), hasDinv_(false), Qpinv_(), hasQpinv_(false), Qvinv_(), hasQvinv_(false)
{;}

LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S, double alp, double beta)
  : S_(S), Ap_(), hasAp_(false), Dinv_(), hasscal_(true), alp_(alp), beta_(beta), hasDinv_(false), Qpinv_(), hasQpinv_(false), Qvinv_(), hasQvinv_(false)
{;}

LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap)
  : S_(S), Ap_(Ap), hasAp_(true), Dinv_(), hasscal_(false), hasDinv_(false), Qpinv_(), hasQpinv_(false), Qvinv_(), hasQvinv_(false)
{;}

LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap, double alp, double beta)
  : S_(S), Ap_(Ap), hasAp_(true), Dinv_(), hasscal_(true), alp_(alp), beta_(beta), hasDinv_(false), Qpinv_(), hasQpinv_(false), Qvinv_(), hasQvinv_(false)
{;}

LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap, double alp, double beta, TSFLinearOperator& Qvinv, TSFLinearOperator& Qpinv)
  : S_(S), Ap_(Ap), hasAp_(true), Dinv_(), hasscal_(true), alp_(alp), beta_(beta), hasDinv_(false), Qpinv_(Qpinv), hasQpinv_(true), Qvinv_(Qvinv), hasQvinv_(true)
{;}

LSCRightOperatorSource::
LSCRightOperatorSource(TSFLinearOperator& S,
                       double alp, double beta, TSFLinearOperator& Qvinv, TSFLinearOperator& Qpinv)
  : S_(S), Ap_(), hasAp_(false), Dinv_(), hasscal_(true), alp_(alp), beta_(beta), hasDinv_(false), Qpinv_(Qpinv), hasQpinv_(true), Qvinv_(Qvinv), hasQvinv_(true)
{;}

TSFReal LSCRightOperatorSource
::getAlpha() const
{
  return alp_;
}

TSFReal LSCRightOperatorSource
::getBeta() const
{
  return beta_;
}

bool LSCRightOperatorSource
::getScal() const
{
  return hasscal_;
}

TSFLinearOperator LSCRightOperatorSource
::getOp() const
{
  return S_;
}

TSFLinearOperator LSCRightOperatorSource
::getAp() const
{

  // if Ap_ already exists (passed in or already built), return it
  if (hasAp_)
    return Ap_;
  
  // if Ap_ doesn't exist yet, need to get it from S
  else
    {
      cerr << "\n LSCRightOperatorSource: in right part of if" << endl;
      
      // Using Ap = C
      Ap_ = S_.getBlock(1,1);
      cerr << "\n LSCRightOperatorSource: got here" << endl;
      
      hasAp_ = true;
      return(Ap_);
      
    }
}

  TSFLinearOperator LSCRightOperatorSource
  ::getDinv() const
  {
//     //SHOULD WE WORRY ABOUT BEING GIVEN A Dinv???
//   // if Dinv_ already exists (passed in or already built), return it

   cerr << "in DiagRightOperatorSource: getDinv()" << endl;
   if (hasDinv_)
     return Dinv_;
//   // if Dinv_ doesn't exist yet, need to get it from F
   else
//     {
      cerr << "DiagRightOperatorSource: in right part of if" << endl;
      TSFLinearOperator C = S_.getBlock(1,1);
      Epetra_CrsMatrix  *C_crs = PetraMatrix::getConcrete(C);

      cerr << "number of diags in F_crs " << C_crs->NumGlobalDiagonals() << endl;

      // get an appropriate vector (with the right map)
      // to hold the extracted diagonals
      // and extract the diagonal from F
      Epetra_Vector Cdiags(C_crs->Map());
      C_crs->ExtractDiagonalCopy(Cdiags);

      // get the reciprocals of the diag values
      Epetra_Vector CdiagsInv(C_crs->Map());
      CdiagsInv.Reciprocal(Cdiags);
      
         // make an epetra matrix for the diagonal matrix
      Epetra_CrsMatrix *Dinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             C_crs->RowMatrixColMap(),
                             C_crs->OperatorDomainMap(),
                             0);

      // Dinv_crs->ReplaceDiagonalValues(FdiagsInv);
      // Note: this is building -Dinv
      int diagSize = CdiagsInv.GlobalLength();
      double value;
      for (int j = 0; j < diagSize; j++)
        {
          value = -1.0*CdiagsInv[j]; // get negative of diags
          Dinv_crs->InsertGlobalValues(j, 1, &value, &j);
        }

      int ierr=Dinv_crs->FillComplete((Epetra_Map) (C_crs->OperatorDomainMap()),
                                          (Epetra_Map) (C_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix FillComplete" << ierr << endl;
        // EPETRA_CHK_ERR(ierr);
      }
      
      // TSF syntax is PetraMatrix(domain, range)
      PetraMatrix* Dinv_petra = new PetraMatrix(C.domain(), C.range());
      Dinv_petra->setPetraMatrix(Dinv_crs,true);
      Dinv_ = Dinv_petra;

      cerr << "number of diags in Dinv_crs " 
           << Dinv_crs->NumGlobalDiagonals() 
           << "\n\n"
           << endl;

      cerr << "\n The matrix is:";
        cerr << Dinv_;
        hasDinv_ = true;
      return(Dinv_);
      
      //   }
  }

  TSFLinearOperator LSCRightOperatorSource
  ::getQvinv() const
  {
//   // if Qvinv_ already exists (passed in or already built), return it

   cerr << "in DiagRightOperatorSource: getQvinv()" << endl;
   if (hasQvinv_)
     return Qvinv_;
//   // if Qvinv_ doesn't exist yet, need to get it from F
   else
//     {
      cerr << "DiagRightOperatorSource: in right part of if" << endl;
      TSFLinearOperator F = S_.getBlock(0,0);
      Epetra_CrsMatrix  *F_crs = PetraMatrix::getConcrete(F);
  
      int numrows = F_crs->NumGlobalDiagonals();
      double N = (double) numrows;

      cerr << "number of diags in F_crs " << numrows << endl;

      // get an appropriate vector (with the right map)
      // to hold the extracted diagonals
      // and extract the diagonal from F
      Epetra_Vector Fdiags(F_crs->Map());
      F_crs->ExtractDiagonalCopy(Fdiags);

         // make an epetra matrix for the diagonal matrix
      Epetra_CrsMatrix *Qvinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             F_crs->RowMatrixColMap(),
                             F_crs->OperatorDomainMap(),
                             0);

      // Note: this is building 1/h^2 I
      int diagSize = Fdiags.GlobalLength();
      double value;
      for (int j = 0; j < diagSize; j++)
        {
          value = 1.0/N; // get negative of diags
          Qvinv_crs->InsertGlobalValues(j, 1, &value, &j);
        }

      int ierr=Qvinv_crs->FillComplete((Epetra_Map) (F_crs->OperatorDomainMap()),
                                          (Epetra_Map) (F_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix FillComplete" << ierr << endl;
        // EPETRA_CHK_ERR(ierr);
      }
      
      Qvinv_crs->OptimizeStorage();

      // TSF syntax is PetraMatrix(domain, range)
      PetraMatrix* Qvinv_petra = new PetraMatrix(F.domain(), F.range());
      Qvinv_petra->setPetraMatrix(Qvinv_crs,true);
      Qvinv_ = Qvinv_petra;

      cerr << "number of diags in Qvinv_crs " 
           << Qvinv_crs->NumGlobalDiagonals() 
           << "\n\n"
           << endl;

      cerr << "\n The matrix is:";
        cerr << Qvinv_;
        hasQvinv_ = true;
      return(Qvinv_);
      
      //   }
  }

  TSFLinearOperator LSCRightOperatorSource
  ::getQpinv() const
  {
//   // if Qpinv_ already exists (passed in or already built), return it

   cerr << "in DiagRightOperatorSource: getQvinv()" << endl;
   if (hasQpinv_)
     return Qpinv_;
//   // if Qpinv_ doesn't exist yet, need to get it from C
   else
//     {
      cerr << "DiagRightOperatorSource: in right part of if" << endl;
      TSFLinearOperator C = S_.getBlock(1,1);
      Epetra_CrsMatrix  *C_crs = PetraMatrix::getConcrete(C);
  
      int numrows = C_crs->NumGlobalDiagonals();
      double N = (double) numrows;

      cerr << "number of diags in C_crs " << numrows << endl;

      // get an appropriate vector (with the right map)
      // to hold the extracted diagonals
      // and extract the diagonal from C
      Epetra_Vector Cdiags(C_crs->Map());
      C_crs->ExtractDiagonalCopy(Cdiags);

         // make an epetra matrix for the diagonal matrix
      Epetra_CrsMatrix *Qpinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             C_crs->RowMatrixColMap(),
                             C_crs->OperatorDomainMap(),
                             0);

      // Qpinv_crs->ReplaceDiagonalValues(CdiagsInv);
      // Note: this is building 1/h^2 I
     
      int diagSize = Cdiags.GlobalLength();
      double value;
      for (int j = 0; j < diagSize; j++)
        {
          value = 1.0/N; // get negative of diags
          Qpinv_crs->InsertGlobalValues(j, 1, &value, &j);
        }

      int ierr=Qpinv_crs->FillComplete((Epetra_Map) (C_crs->OperatorDomainMap()),
                                          (Epetra_Map) (C_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix FillComplete" << ierr << endl;
        // EPETRA_CHK_ERR(ierr);
      }
      
      Qpinv_crs->OptimizeStorage();

      // TSF syntax is PetraMatrix(domain, range)
      PetraMatrix* Qpinv_petra = new PetraMatrix(C.domain(), C.range());
      Qpinv_petra->setPetraMatrix(Qpinv_crs,true);
      Qpinv_ = Qpinv_petra;

      cerr << "number of diags in Qpinv_crs " 
           << Qpinv_crs->NumGlobalDiagonals() 
           << "\n\n"
           << endl;

      cerr << "\n The matrix is:";
        cerr << Qpinv_;
        hasQpinv_ = true;
      return(Qpinv_);
      
      //   }
  }

string LSCRightOperatorSource::toString() const 
{
	return "LSC operator source";
}
