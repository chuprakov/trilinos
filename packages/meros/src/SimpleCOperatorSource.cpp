   #include "SimpleCOperatorSource.h"
   #include "RightBlockNSOperatorSource.h"
   #include "TSFOperatorSourceBase.h"
   #include "TSFBlockLinearOperator.h"
   #include "Epetra_CrsMatrix.h"
   #include "Epetra_Vector.h"
   #include "PetraMatrix.h"

  using namespace TSF;
  using namespace SPP;
  using std::string;

  SimpleCOperatorSource::SimpleCOperatorSource(TSFLinearOperator& S)
  : S_(S), Dinv_(), hasDinv_(false)
  {;}

  TSFLinearOperator SimpleCOperatorSource
  ::getOp() const
 {
  return S_;
 }

  TSFLinearOperator SimpleCOperatorSource
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

       cerr << "\n Simple C " << endl;

  TSFVector vec = F.range().createMember();
  TSFVector posvec = F.range().createMember();
  TSFVector finalvec = F.range().createMember();
  int numRows = F.range().dim();
  int len;

  for (int i = 0; i<numRows; i++)
    {
    TSFArray<int> ind;
    TSFArray<TSFReal> aij;
    F.getRow(i, ind, aij);
    len = aij.length();
    for (int j = 0; j<len; j++)
      {
	vec[j] = aij[j];
      }
    //posvec = vec.abs();
    TSFReal sum  = vec.abs().sumElements();
    finalvec[i] = sum;
    }
  cerr << "\nfinalvec is:" << finalvec;
   Epetra_CrsMatrix *Dinv_crs = 
        new Epetra_CrsMatrix(Copy, 
                             F_crs->RowMatrixColMap(),
                             F_crs->OperatorDomainMap(),
                             0);

 for (int j = 0; j < numRows; j++)
        {
          double value = -1.0/finalvec[j]; // get negative of diags
          Dinv_crs->InsertGlobalValues(j, 1, &value, &j);
        }

      int ierr=Dinv_crs->TransformToLocal((Epetra_Map *) &(F_crs->OperatorDomainMap()),
                                          (Epetra_Map *) &(F_crs->RowMatrixColMap()));
      if (ierr!=0) {
        cerr <<"Error in Epetra_CrsMatrix TransformToLocal" << ierr << endl;
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

      cerr << "\n The matrix is:";
        cerr << Dinv_;
      hasDinv_ = true;
      return(Dinv_);
      
    }
  }

  string SimpleCOperatorSource::toString() const 
  {
	return "Simple C operator source";
  }
