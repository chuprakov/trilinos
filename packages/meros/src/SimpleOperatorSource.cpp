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

   #include "SimpleOperatorSource.h"
   #include "RightBlockNSOperatorSource.h"
   #include "TSFOperatorSourceBase.h"
   #include "TSFBlockLinearOperator.h"
   #include "Epetra_CrsMatrix.h"
   #include "Epetra_Vector.h"
   #include "PetraMatrix.h"

  using namespace TSF;
  using namespace Meros;
  using std::string;

  SimpleOperatorSource::SimpleOperatorSource(TSFLinearOperator& S)
  : S_(S), Dinv_(), hasDinv_(false)
  {;}

  TSFLinearOperator SimpleOperatorSource
  ::getOp() const
 {
  return S_;
 }

  TSFLinearOperator SimpleOperatorSource
  ::getDinv() const
  {
  // if Dinv_ already exists (passed in or already built), return it

  cerr << "in DiagRightOperatorSource: getDinv()" << endl;
  if (hasDinv_)
    return Dinv_;
  // if Dinv_ doesn't exist yet, need to get it from F
  else
    {
      cerr << "DiagRightOperatorSource: in right part of if" << endl;
      TSFLinearOperator F = S_.getBlock(0,0);
      Epetra_CrsMatrix  *F_crs = PetraMatrix::getConcrete(F);

      cerr << "number of diags in F_crs " << F_crs->NumGlobalDiagonals() << endl;

      // get an appropriate vector (with the right map)
      // to hold the extracted diagonals
      // and extract the diagonal from F
      Epetra_Vector Fdiags(F_crs->Map());
      F_crs->ExtractDiagonalCopy(Fdiags);

      // get the reciprocals of the diag values
      Epetra_Vector FdiagsInv(F_crs->Map());
      FdiagsInv.Reciprocal(Fdiags);
      
      // make an epetra matrix for the diagonal matrix
      Epetra_CrsMatrix *Dinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             F_crs->RowMatrixColMap(),
                             F_crs->OperatorDomainMap(),
                             0);

      // Dinv_crs->ReplaceDiagonalValues(FdiagsInv);
      // Note: this is building -Dinv
      int diagSize = FdiagsInv.MyLength();
      cerr << "\n diagsize is: " << diagSize;
      double value;
      for (int j = 0; j < diagSize; j++)
        {
          value = -1.0*FdiagsInv[j]; // get negative of diags
          Dinv_crs->InsertMyValues(j, 1, &value, &j);
        }

      int ierr=Dinv_crs->FillComplete((Epetra_Map) (F_crs->OperatorDomainMap()),
                                          (Epetra_Map) (F_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix FillComplete" << ierr << endl;
        // EPETRA_CHK_ERR(ierr);
      }
      
      // TSF syntax is PetraMatrix(domain, range)
      PetraMatrix* Dinv_petra = new PetraMatrix(F.domain(), F.range());
      Dinv_petra->setPetraMatrix(Dinv_crs,true);
      Dinv_ = Dinv_petra;

      cerr << "number of diags in Dinv_crs " 
           << Dinv_crs->NumGlobalDiagonals() 
           << "\n\n"
           << endl;

      //      cerr << "\n The matrix is:";
      // cerr << Dinv_;
      hasDinv_ = true;
      return(Dinv_);
      
    }
  }

  string SimpleOperatorSource::toString() const 
  {
	return "Simple operator source";
  }
