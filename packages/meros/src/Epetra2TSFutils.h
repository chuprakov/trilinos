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


