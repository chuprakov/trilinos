#include "DiagRightOperatorSource.h"
#include "RightBlockNSOperatorSource.h"
#include "TSFOperatorSourceBase.h"
#include "TSFBlockLinearOperator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "PetraMatrix.h"

using namespace TSF;
using namespace SPP;
using std::string;


DiagRightOperatorSource::
DiagRightOperatorSource(TSFLinearOperator& S)
 	: S_(S), Dinv_(), hasDinv_(false)
{;}

DiagRightOperatorSource::
DiagRightOperatorSource(TSFLinearOperator& S,
                        TSFLinearOperator& D)
 	: S_(S), Dinv_(D), hasDinv_(true)
{;}


TSFLinearOperator DiagRightOperatorSource
::getOp() const
{
  return S_;
}

TSFLinearOperator DiagRightOperatorSource
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
      
      // cerr << "F diags" << endl;
      // cerr << FdiagsInv << endl;
      // cerr << Fdiags << endl;

      // make an epetra matrix for the diagonal matrix
      Epetra_CrsMatrix *Dinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             F_crs->RowMatrixColMap(),
                             F_crs->OperatorDomainMap(),
                             0);

      // Dinv_crs->ReplaceDiagonalValues(FdiagsInv);
      // Note: this is building -Dinv
      int diagSize = FdiagsInv.GlobalLength();
      double value;
      for (int j = 0; j < diagSize; j++)
        {
          value = -1.0*FdiagsInv[j]; // get negative of diags
          Dinv_crs->InsertGlobalValues(j, 1, &value, &j);
        }
      int ierr=Dinv_crs->TransformToLocal((Epetra_Map *) &(F_crs->OperatorDomainMap()),
                                          (Epetra_Map *) &(F_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix TransformToLocal" << ierr << endl;
        // EPETRA_CHK_ERR(ierr);
      }
      
      PetraMatrix* Dinv_petra = new PetraMatrix(F.range(), F.domain());
      Dinv_petra->setPetraMatrix(Dinv_crs,true);
      Dinv_ = Dinv_petra;

      cerr << "number of diags in Dinv_crs " 
           << Dinv_crs->NumGlobalDiagonals() 
           << "\n\n"
           << endl;

      hasDinv_ = true;
      return(Dinv_);
      
    }
}


string DiagRightOperatorSource::toString() const 
{
	return "Diag operator source";
}
