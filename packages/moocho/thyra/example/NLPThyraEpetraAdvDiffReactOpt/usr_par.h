#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"

#include <iostream>

/** \class Usr_Par
    \brief Computes and stores several quantities used in the local FE assembly. There is only
    one method \e Print, in addition to the constructor. All computations are performed
    during construction.
*/
            
class Usr_Par {
public:

  Usr_Par();

  Epetra_SerialDenseMatrix Nodes;
  Epetra_SerialDenseVector Weights;

  Epetra_SerialDenseMatrix N;

  Epetra_SerialDenseMatrix Nx1;

  Epetra_SerialDenseMatrix Nx2;

  Epetra_SerialDenseMatrix S1;
  Epetra_SerialDenseMatrix S2;
  Epetra_SerialDenseMatrix S3;

  Epetra_SerialDenseVector Nw;

  Epetra_SerialDenseMatrix NNw;

  Epetra_SerialDenseMatrix * NNNw;

  Epetra_SerialDenseMatrix * NdNdx1Nw;

  Epetra_SerialDenseMatrix * NdNdx2Nw;

  ~Usr_Par() {
    delete [] NNNw;
    delete [] NdNdx1Nw;
    delete [] NdNdx2Nw;
  }
  
  void Print(ostream& os) const;
};
