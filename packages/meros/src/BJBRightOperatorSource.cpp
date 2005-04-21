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

#include "BJBRightOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFBlockLinearOperator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "PetraMatrix.h"

using namespace TSF;
using namespace Meros;
using std::string;


BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S)
  : S_(S), Ap_(), hasAp_(false), Dinv_()
{;}

BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S, double alp, double beta)
  : S_(S), Ap_(), hasAp_(false), Dinv_(), hasscal_(true), alp_(alp), beta_(beta)
{;}

BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap)
  : S_(S), Ap_(Ap), hasAp_(true), Dinv_(), hasscal_(false)
{;}

BJBRightOperatorSource::
BJBRightOperatorSource(TSFLinearOperator& S,
                       TSFLinearOperator& Ap, double alp, double beta)
  : S_(S), Ap_(Ap), hasAp_(true), Dinv_(), hasscal_(true), alp_(alp), beta_(beta)
{;}

TSFReal BJBRightOperatorSource
::getAlpha() const
{
  return alp_;
}

TSFReal BJBRightOperatorSource
::getBeta() const
{
  return beta_;
}

bool BJBRightOperatorSource
::getScal() const
{
  return hasscal_;
}

TSFLinearOperator BJBRightOperatorSource
::getOp() const
{
  return S_;
}

TSFLinearOperator BJBRightOperatorSource
::getAp() const
{

  // if Ap_ already exists (passed in or already built), return it
  if (hasAp_)
    return Ap_;
  
  // if Ap_ doesn't exist yet, need to get it from S
  else
    {
      cerr << "BJBRightOperatorSource: in right part of if" << endl;
      
      // Using Ap = C
      Ap_ = S_.getBlock(1,1);
      cerr << "BJBRightOperatorSource: got here" << endl;
      
      hasAp_ = true;
      return(Ap_);
      
    }
}

  TSFLinearOperator BJBRightOperatorSource
  ::getDinv() const
  {
//     //SHOULD WE WORRY ABOUT BEING GIVEN A Dinv???
//   // if Dinv_ already exists (passed in or already built), return it

//   cerr << "in DiagRightOperatorSource: getDinv()" << endl;
//   if (hasDinv_)
//     return Dinv_;
//   // if Dinv_ doesn't exist yet, need to get it from F
//   else
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
	//      hasDinv_ = true;
      return(Dinv_);
      
      //   }
  }

string BJBRightOperatorSource::toString() const 
{
	return "BJB operator source";
}
